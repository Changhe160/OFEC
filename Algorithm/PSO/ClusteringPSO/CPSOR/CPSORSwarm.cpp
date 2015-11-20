#include "CPSORSwarm.h"

#ifdef OFEC_DEMON
#include "../../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

CPSORSwarm::CPSORSwarm(ParamMap &v):Swarm(),MultiPopulationCont<CPSORSubSwarm>((v[param_subPopSize]),(v[param_overlapDgre])),\
	m_subSize(7){
	m_name="ALG_CPSOR";
	if(IS_PROBLEM_NAME(Global::ms_curProId,"DYN_CONT_MovingPeak")){
		m_diversityDegree=1-exp(-0.2*pow(CAST_PROBLEM_DYN_CONT->getNumberofPeak(),0.45));
		m_initialSize =300*(1-exp(-0.33*pow(CAST_PROBLEM_DYN_CONT->getNumberofPeak(),0.5)));
	}else{
		m_diversityDegree=0.8;
		m_initialSize =100;
	}
	initialize();
	m_algPar<<"Initial population size: "<<m_initialSize<<"; sub-population size: "<<m_subSize<<"; Overlap degree:"<<v[param_overlapDgre]<<";"<<"; Diversity degree:"<<m_diversityDegree<<";";
    
}
void CPSORSwarm::initialize(){
	
	int size=0;
	for(size_t i=0;i<m_subPop.size();i++)
		size+=m_subPop[i]->getPopSize();
	m_popsize=m_initialSize-size; //the number of particles to be added
	m_pop.clear();
	for(int i=0;i<m_popsize;i++)
		m_pop.push_back(move(unique_ptr<CPSORParticle>(new CPSORParticle())));
	if(CAST_PROBLEM_DYN->predictChange(m_popsize)){
		Population::initialize(false,false,true);
	}else{
		Population::initialize(false,true,true);
	}
}

void CPSORSwarm::createSubswarms( ){

	m_clst.setNormalizationFlag(false);
	m_clst.setSpace(0);
	m_clst.initialize(m_pop,m_popsize);
	m_clst.roughClustering(m_subSize);

	for(int k=0;k< m_clst.getSize();k++){			
		CPSORSubSwarm *s=new CPSORSubSwarm(m_clst[k]);
		s->updateCurRadius(true);
		addPopulation(*s);
	}
	m_clst.clear_();
	remove(m_popsize);

	measureMultiPop();
}

ReturnFlag CPSORSwarm::run_(){
	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	createSubswarms();
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
		//cout<<Global::msp_global->mp_problem->getEvaluations()<<" "<<getNumPops()<<" "<<m_subPop[findBestPop(1)]->m_best[0]->obj(0)<<endl;
		measureMultiPop();


		if(m_subPop.size()>1){
			while(removeOverlapping()!=-1);
			for(auto it=m_subPop.begin();it!=m_subPop.end();++it) (*it)->checkOverCrowd(m_subSize);
				
			for(decltype(m_subPop.size()) i=0;i<m_subPop.size();i++){
				if(m_subPop[i]->isConverged(0.0001)){
					deletePopulation(i);
					i--;
				}
			}
		}		

		if(Global::msp_global->m_totalNumIndis<m_initialSize*m_diversityDegree){
			initialize();
			if(ifTerminating()) break;
			createSubswarms();
		}

	}
	return Return_Terminate;
}