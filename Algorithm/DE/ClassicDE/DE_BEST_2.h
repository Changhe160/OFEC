#ifndef DE_BEST_2_H
#define DE_BEST_2_H

#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

#include "../DEIndividual.h"
#include "../DEPopulation.h"
class DEBest2:public DEPopulation<CodeVReal,DEIndividual>{
public:
	DEBest2(ParamMap &v):DEPopulation(v[param_popSize]){
		setMutationStrategy(DE_best_2);
	}
	ReturnFlag run_(){
		ReturnFlag rf=Return_Normal;

		#ifdef OFEC_CONSOLE
		if(Global::msp_global->m_runId==0){
			mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
			mMultiModal::getPopInfor()->setOutProgFlag(true);
		}
		#endif // OFEC_CONSOLE

		while(!ifTerminating()){
			/*g_mutexStream.lock();
			cout<<Global::msp_global->m_runId<<" "<<Global::msp_global->mp_problem->getEvaluations()<<" "<<m_best[0]->obj(0)<<endl;
			g_mutexStream.unlock();*/
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
	 
		 return rf;
	}
};
#endif