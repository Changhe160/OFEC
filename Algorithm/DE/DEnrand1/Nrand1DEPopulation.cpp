#include "Nrand1DEPopulation.h"
#include "../../../Measure/mMultiObj.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif
Nrand1DEPopulation::Nrand1DEPopulation(ParamMap &v) :DEPopulation<CodeVReal,Nrand1DEIndividual>(v[param_popSize],true)
{
	m_CR=0.9;
	m_F=0.5;
}

void Nrand1DEPopulation::mutate(const int idx)
{
	int nearest=findNearest(idx,0,0,4);
	vector<int>a(this->m_popsize);
	Global::msp_global->initializeRandomArray<vector<int>>(a,this->m_popsize);
    int j=0;
    while(a[j]!=idx){j++;}
    int r1,r2;
    r1=a[(j+1)%this->m_popsize];
    r2=a[(j+2)%this->m_popsize];
	this->m_pop[idx]->mutate(m_F,&this->m_pop[nearest]->self(),&this->m_pop[r1]->self(),&this->m_pop[r2]->self());
}


ReturnFlag Nrand1DEPopulation::run_()
{
	ReturnFlag rf=Return_Normal;

	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	while(!ifTerminating())
	{
		#ifdef OFEC_DEMON
			vector<Algorithm*> vp;	
			vp.push_back(this);	
			msp_buffer->updateBuffer_(&vp);
		#endif
			
		rf=evolve();
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

bool Nrand1DEPopulation::ifTerminating(){
	
	#ifdef OFEC_DEMON
		return Algorithm::ifTerminating();
	#endif

	#ifdef OFEC_CONSOLE
	
	if(Global::msp_global->mp_problem->m_name.find("FUN_")!=string::npos){
		if(CAST_PROBLEM_CONT->getGOpt().getNumGOptFound()==CAST_PROBLEM_CONT->getGOpt().getNumOpt()||Algorithm::ifTerminating()){
				
			if(mMultiModal::getPopInfor()){
				int peaksf;		
				peaksf=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
				mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
					Global::msp_global->m_totalNumIndis,1,peaksf,\
					CAST_PROBLEM_CONT->getGOpt().getNumOpt(),0,0,0,0,\
					0,0,CAST_PROBLEM_CONT->getGOpt().isAllFound());
			}
				
			return true;
		}
	}
	return false;
	#endif
}