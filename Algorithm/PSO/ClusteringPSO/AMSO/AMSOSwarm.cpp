#include "AMSOSwarm.h"

#ifdef OFEC_DEMON
#include "../../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

AMSOSwarm::AMSOSwarm(ParamMap &v):Swarm<CodeVReal,AMSOParticle>(),MultiPopulationCont<AMSOSubSwarm>((v[param_subPopSize]),(v[param_overlapDgre])),\
	m_subSize(7),m_maxIndis(300),m_minIndis(70),m_stepIndis(10),m_gap(1500),m_alpha(-0.002),m_offPeak(3),\
	m_convergedPops(0),m_preIndis(0),m_nextIndis(0),m_bufferTimes(0),m_maxPops(0)
{
	m_name="ALG_AMSO";
	for(int i=0;i<this->m_popsize;i++)
		this->m_pop[i]->setVmax(-this->m_initialRadius,this->m_initialRadius);
	m_initialSize=100;
	initialize();
	m_algPar<<"Initial population size: "<<m_initialSize<<"; sub-population size: "<<m_subSize<<"; Overlap degree:"<<v[param_overlapDgre]<<";";
}

AMSOSwarm::~AMSOSwarm()
{
	mv_converge.clear();
	while(!mq_pops.empty())
		 mq_pops.pop();
}

void AMSOSwarm::initialize()
{
	int size=0;
	for(size_t i=0;i<m_subPop.size();i++)
		size+=m_subPop[i]->getPopSize();
	size+=mv_converge.size();
	m_popsize=m_initialSize-size; //the number of particles to be added
	m_pop.clear();
	for(int i=0;i<m_popsize;i++)
		m_pop.push_back(move(unique_ptr<AMSOParticle>(new AMSOParticle())));
	if(CAST_PROBLEM_DYN->predictChange(m_popsize)){
		Population<CodeVReal,AMSOParticle>::initialize(false,false,true);
	}else{
		Population<CodeVReal,AMSOParticle>::initialize(false,true,true);
	}
	
	for(auto &i:mv_converge){
		if(i->getNumDim()>GET_NUM_DIM)	i->decreaseDimension();
		if(i->getNumDim()<GET_NUM_DIM) i->increaseDimension();
		i->updateMemory();
	}
	PopulationCont::add(mv_converge);
	mv_converge.clear();
}

void AMSOSwarm::createSubswarms()
{
	m_clst.setNormalizationFlag(false);
	m_clst.setSpace(0);
	m_clst.initialize(m_pop,m_popsize);
	m_clst.roughClustering(m_subSize);

	for(int k=0;k<m_clst.getSize();k++){			
		AMSOSubSwarm *s=new AMSOSubSwarm(m_clst[k]);
		s->updateCurRadius(true);
		addPopulation(*s);
	}
	m_clst.clear_();
	remove(m_popsize);

	measureMultiPop();
	PopInfor infor(Global::msp_global->mp_problem->getEvaluations(),m_subPop.size()+m_convergedPops);
    if(mq_pops.size()>0&&mq_pops.back().m_evals-mq_pops.front().m_evals>=m_gap){
		// set mq_pops to be initial status if the difference between the fitness evaluations of the back and the front member >=mc_gap
         mq_pops.pop();
     }
     //pop in the populations infor to mp_pops
     mq_pops.push(infor);
}

void AMSOSwarm::getNextIndis(PopInfor &infor)
{
	if(m_subPop.size()==0||mq_pops.size()>0&&infor.m_evals-mq_pops.front().m_evals>=m_gap)
	{
		double ratio=(mq_pops.back().m_pops-mq_pops.front().m_pops)/(double)(mq_pops.back().m_evals-mq_pops.front().m_evals);
		// calculate the pop number decreasing ratio
		if(m_subPop.size()==0||ratio>m_alpha)
		{
			if(m_subPop.size()==0||m_bufferTimes<=1)
			{
				//the number of individuals changed for the last check point
				if(m_subPop.size()==0) m_nextIndis=m_minIndis;
				else m_nextIndis=m_preIndis;
			}
			else
			{
				if(infor.m_pops-m_maxPops>0)
					m_nextIndis=m_preIndis+m_stepIndis*(infor.m_pops-m_maxPops);
				else if(infor.m_pops-m_maxPops<-m_offPeak)
					m_nextIndis=m_preIndis-m_stepIndis*(m_maxPops-infor.m_pops);
				else
					m_nextIndis=m_preIndis;
			}

			increaseDiversity();
              
			if(m_nextIndis!=m_preIndis){
				m_bufferTimes=1;
				m_maxPops=infor.m_pops;
			}else {
				m_bufferTimes++;
				if(m_maxPops<infor.m_pops) m_maxPops=infor.m_pops;
			}

			m_preIndis=m_initialSize;

			infor.m_evals=Global::msp_global->mp_problem->getEvaluations();
			infor.m_pops=m_subPop.size(); 

			while(!mq_pops.empty()){
				mq_pops.pop();
			}
			m_convergedPops=0;

			if(mv_converge.size()>0) mv_converge.clear();
		}
	}
}

void AMSOSwarm::increaseDiversity()
{
	if(m_nextIndis>m_maxIndis) m_nextIndis=m_maxIndis;
	if(m_nextIndis<m_minIndis) m_nextIndis=m_minIndis;
	int number=0;
	for(size_t i=0;i<m_subPop.size();i++)
		number+=m_subPop[i]->getPopSize();
	if(m_nextIndis<number)   return;
	number+=mv_converge.size();
	number=m_nextIndis-number;
	if(number>0)
	{
		m_initialSize=m_nextIndis;
		initialize();
		createSubswarms();
	}
}

ReturnFlag AMSOSwarm::run_(){
	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	createSubswarms();
	m_preIndis=m_initialSize;
	ReturnFlag r_flag=Return_Normal;
	while(!ifTerminating()){
		for(auto &swarm:m_subPop){
			if(swarm->m_popsize==0) continue;
			r_flag=swarm->evolve();
			if(r_flag==Return_Terminate) break;

			handleReturnFlagAll(r_flag);
			HANDLE_RETURN_FLAG(r_flag)
			
			#ifdef OFEC_DEMON
					vector<Algorithm*> vp;
					for(auto &it:m_subPop){
						vp.push_back(it.get());
					}
					msp_buffer->updateBuffer_(&vp);
			#endif
		}
		if(r_flag==Return_Terminate) break;

		if(m_subPop.size()>1)
		{				
			while(removeOverlapping()!=-1);
			for(auto it=m_subPop.begin()+1;it!=m_subPop.end();++it) (*it)->checkOverCrowd(m_subSize);
				
			for(decltype(m_subPop.size()) i=1;i<m_subPop.size();i++)
			{
				if(m_subPop[i]->isConverged(0.0001))
				{
					for(auto &i:mv_converge)
					{
						if(i->getNumDim()>GET_NUM_DIM)	i->decreaseDimension();
						if(i->getNumDim()<GET_NUM_DIM) i->increaseDimension();
						i->updateMemory();
					}
					mv_converge.push_back(move(unique_ptr<AMSOParticle>(new AMSOParticle(*(m_subPop[i]->m_best[0])))));
					deletePopulation(i);
					i--;
					m_convergedPops++;
				}
			}				
		}
		measureMultiPop();

		//create infor to record current status
		PopInfor infor(Global::msp_global->mp_problem->getEvaluations(),m_subPop.size()+m_convergedPops);

		getNextIndis(infor);

		if(mq_pops.size()>0&&mq_pops.back().m_evals-mq_pops.front().m_evals>=m_gap){
		// set mq_pops to be initial status if the difference between the fitness evaluations of the back and the front member >=mc_gap
			mq_pops.pop();
		}
		mq_pops.push(infor);
	}
	return Return_Terminate;
}