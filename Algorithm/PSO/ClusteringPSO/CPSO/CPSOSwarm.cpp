#include "CPSOSwarm.h"

#ifdef OFEC_DEMON
#include "../../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

CPSOSwarm::CPSOSwarm(ParamMap &v):Algorithm(-1,"ALG_CPSO"),MultiPopulationCont<CPSOSubSwarm>((v[param_subPopSize]),(v[param_overlapDgre])),\
	m_initialSize(v[param_popSize]),m_subSize(v[param_subPopSize]){
	
}
CPSOSwarm::~CPSOSwarm(){
	mv_converge.clear();
};
void CPSOSwarm::initialize(){

	int size=m_initialSize-mv_converge.size()-m_subPop.size()+1;
	CPSOSubSwarm *s=nullptr;
	if(CAST_PROBLEM_DYN->predictChange(size))	s=new CPSOSubSwarm(size,false);
	else	s=new CPSOSubSwarm(size,true);

	for(auto &i:mv_converge){
		if(i->getNumDim()>GET_NUM_DIM)	i->decreaseDimension();
		if(i->getNumDim()<GET_NUM_DIM) i->increaseDimension();
		i->updateMemory();
	}
	s->add(mv_converge);

	for(auto &swarm:m_subPop){
		if(swarm->m_popsize==0) continue;
		s->add(new CPSOParticle(*(swarm->m_best[0])));
	}
	m_subPop.clear();
	mv_converge.clear();
	addPopulation(*s);
	
}

void CPSOSwarm::createSubswarms( ){

	m_clst.setNormalizationFlag(false);
	m_clst.setSpace(0);
	m_clst.initialize(m_subPop[0]->m_pop,m_subPop[0]->getPopSize());
	m_clst.roughClustering(m_subSize);

	for(int k=0;k< m_clst.getSize();k++){			
		CPSOSubSwarm *s=new CPSOSubSwarm(m_clst[k]);
		s->updateCurRadius(true);
		addPopulation(*s);
		
	}
	m_clst.clear_();
	m_subPop[0]->remove(m_subPop[0]->getPopSize());

	measureMultiPop();
}

ReturnFlag CPSOSwarm::run_(){
	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	initialize();
	createSubswarms();
	ReturnFlag r_flag=Return_Normal;
	while(!ifTerminating()){
		bool flag=false;
		for(auto &swarm:m_subPop){
			if(swarm->m_popsize==0) continue;
			flag=false;
			Solution<CodeVReal> x(*swarm->m_best[0]);
			vector<double> objOld=x.obj(); 
			r_flag=swarm->evolve();
			if(r_flag==Return_Terminate) break;
			handleReturnFlagAll(r_flag);
			HANDLE_RETURN_FLAG(r_flag)

			if(!CAST_PROBLEM_DYN->predictChange(1)){
				r_flag=x.evaluate(true);
			}
			if(x.obj(0)!=objOld[0])	flag=true;
			
			handleReturnFlagAll(r_flag);
			HANDLE_RETURN_FLAG(r_flag)
			
			#ifdef OFEC_DEMON
					vector<Algorithm*> vp;
					for(auto &it:m_subPop){
						vp.push_back(it.get());
					}
					msp_buffer->updateBuffer_(&vp);
			#endif
			if(flag) break;
		}
		if(r_flag==Return_Terminate) break;
		//cout<<Global::msp_global->mp_problem->getEvaluations()<<" "<<getNumPops()<<" "<<m_subPop[findBestPop(1)]->m_best[0]->obj(0)<<endl;
		measureMultiPop();

		if(flag){
			initialize();
			if(ifTerminating()) break;
			createSubswarms();
		}else{
			if(Global::msp_global->m_totalNumIndis>=m_subSize){				
				while(removeOverlapping()!=-1);
				for(auto it=m_subPop.begin()+1;it!=m_subPop.end();++it) (*it)->checkOverCrowd(m_subSize);
				
				for(decltype(m_subPop.size()) i=1;i<m_subPop.size();i++){
					if(m_subPop[i]->isConverged(0.0001)){
						if(m_initialSize-mv_converge.size()-m_subPop.size()+1>1){
							mv_converge.push_back(move(unique_ptr<CPSOParticle>(new CPSOParticle(*(m_subPop[i]->m_best[0])))));
						}
						deletePopulation(i);
						i--;
					}
				}				
			}else{
				int num=m_subSize-Global::msp_global->m_totalNumIndis;
				if(CAST_PROBLEM_DYN->predictChange(num))			m_subPop[0]->add(num,false,false);
				else m_subPop[0]->add(num,true,false);
			}

		}

	}
	return Return_Terminate;
}

int CPSOSwarm::removeOverlapping(){
	for(unsigned i=1;i<this->m_subPop.size();i++){
		if(this->m_subPop[i]->m_popsize==0) continue;
		for(unsigned j=i+1;j<this->m_subPop.size();j++){	
			if(this->m_subPop[j]->m_popsize==0) continue;
			double dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_center);
			if(dist<this->m_subPop[i]->m_initialRadius||dist<this->m_subPop[j]->m_initialRadius){
				int c1=0,c2=0;
				for(int k=0;k<this->m_subPop[j]->m_popsize;k++){
					dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c1++;
				}
				for(int k=0;k<this->m_subPop[i]->m_popsize;k++){
					dist=this->m_subPop[j]->m_center.getDistance(this->m_subPop[i]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c2++;
				}
				if(c1>this->m_subPop[j]->m_popsize*m_overlapDegree&&c2>this->m_subPop[i]->m_popsize*m_overlapDegree){
					int idx=-1;
					if(*this->m_subPop[i]>(*this->m_subPop[j])){	
						this->m_subPop[i]->add(*this->m_subPop[j]);
						this->deletePopulation(j);
						idx=j;
					}else if(*this->m_subPop[i]<(*this->m_subPop[j])){
						this->m_subPop[j]->add(*this->m_subPop[i]);
						this->deletePopulation(i);
						idx=i;
					}else{
						throw myException("TODO: comparison between two populations@removeOverlapping()");
					}
					return idx;
				}
			}
		}
	}
	return -1;
}