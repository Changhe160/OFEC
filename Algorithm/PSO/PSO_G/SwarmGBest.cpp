#include "SwarmGBest.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

SwarmGBest::SwarmGBest() :Swarm<CodeVReal,Particle>()
{
}

SwarmGBest::SwarmGBest(ParamMap &v) :Swarm<CodeVReal,Particle>(v[param_popSize],true)
{
	
}

SwarmGBest::SwarmGBest(int popsize, bool mode):Swarm<CodeVReal,Particle>(popsize,mode)
{
}

SwarmGBest::SwarmGBest(const Solution<CodeVReal> & center, double radius, int rPopsize,bool mode) :Swarm<CodeVReal,Particle>(center,radius,rPopsize,mode)
{
}


ReturnFlag SwarmGBest::evolve(){
	if(this->m_popsize<1) return Return_Normal;
	ReturnFlag r_flag=Return_Normal;

	for(int i=0;i<this->m_popsize;i++){				  
		r_flag=this->m_pop[i]->move(neighborBest(i),m_W,m_C1,m_C2);//			
		if(this->m_pop[i]->self()>this->m_pop[i]->representative()){
			this->m_pop[i]->representative()=this->m_pop[i]->self();
			this->updateBestArchive(this->m_pop[i]->self());
		}
		if(r_flag!=Return_Normal) break;
	}
	if(r_flag==Return_Normal){
		this->m_iter++;
	}
	return r_flag;
}

ReturnFlag SwarmGBest::run_()
{
	ReturnFlag rf=Return_Normal;

	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
	if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	int maxG = Global::g_arg[param_maxEvals] / m_popsize;

	while(!ifTerminating()){
		g_mutexStream.lock();
		//cout<<Global::msp_global->m_runId<<" "<<Global::msp_global->mp_problem->getEvaluations()<<" "<<m_best[0]->obj(0)<<endl;
		g_mutexStream.unlock();
		rf=this->evolve();

		#ifdef OFEC_DEMON
			vector<Algorithm*> vp;	
			vp.push_back(this);	
			msp_buffer->updateBuffer_(&vp);
		#endif
		if(rf!=Return_Normal) handleReturnFlag(rf);

		#ifdef OFEC_CONSOLE
		if(mMultiModal::getPopInfor()){
			int peaksf;		
			peaksf=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
			mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
				Global::msp_global->m_totalNumIndis,1,peaksf,\
				CAST_PROBLEM_CONT->getGOpt().getNumOpt(),0,0,0,0,\
				0,0,CAST_PROBLEM_CONT->getGOpt().isAllFound());
		}
		#endif

		if(rf==Return_Terminate) break;
    }
	cout << "run" << Global::msp_global->m_runId << endl;
     return rf;
}