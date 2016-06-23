#include "CRDEPopulation.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

CRDEPopulation::CRDEPopulation(ParamMap &v) :DEPopulation<CodeVReal,DEIndividual>(v[param_popSize],true)
{
	m_F=0.5;
	m_CR=0.9;
	m_mutStrategy=DE_rand_1;
}


ReturnFlag CRDEPopulation::evolve()
{
	if(this->m_popsize<4){
       throw myException("the population size cannot be smaller than 5@DEPopulation<TypeDEIndi>::evolve()");
	}

	ReturnFlag rf=Return_Normal;

    for(int i=0;i<this->m_popsize;i++){
        mutate(i);
        this->m_pop[i]->recombine(m_CR);
  
		rf=this->m_pop[i]->m_pu.evaluate();
		if(rf!=Return_Normal) return rf;

        int idx=this->findNearest(i,0,0,1);
	
		if(this->m_pop[i]->m_pu>this->m_pop[idx]->self()){
			this->m_pop[idx]->self()=this->m_pop[i]->m_pu;
		}
		this->updateBestArchive(DEIndividual(this->m_pop[i]->m_pu));
	}
    this->m_iter++;
	return rf;
}

ReturnFlag CRDEPopulation::run_()
{
	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	ReturnFlag rf=Return_Normal;
	while(!ifTerminating())
	{
	//	cout<<Global::msp_global->mp_problem->getEvaluations()<<endl;

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

bool CRDEPopulation::ifTerminating(){
	
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