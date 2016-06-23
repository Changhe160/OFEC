#include "CMAES_.h"
#include "cmaes.h"
#include "cmaes_interface.h"
#include "../../../Global/global.h"
#include "../../../Measure/mMultiModal.h"
#include "../../../Measure/mSingleObj.h"
#include "../../../Problem/FunctionOpt/FFreePeak_OnePeak.h"

CMAES::CMAES(ParamMap &v):Algorithm(-1, "CMAES"), m_pop(v[param_popSize])
{

	 initialsFilePathName =(string) v[param_workingDir];
	 signalsFilePathName= (string) v[param_workingDir];

	 initialsFilePathName+= "Algorithm/CMAES/initials.par";
	 signalsFilePathName+= "Algorithm/CMAES/signals.par";
}

CMAES::~CMAES()
{
     //dtor
	
}

void CMAES::copy(double*x, CodeVReal &vx) {
	int numDim = Global::msp_global->mp_problem->getNumDim();
	
	for (int i = 0; i<numDim; ++i) {
		vx[i] = x[i];
	}
}
double CMAES::fitCompute(Solution<CodeVReal> &s)
{	
	s.evaluate(true);
	if(Global::msp_global->mp_problem->getOptType()==MIN_OPT)    return s.data().m_obj[0];
	else return -s.data().m_obj[0];
}

ReturnFlag CMAES::run_(){
#ifdef OFEC_CONSOLE
	if (Global::msp_global->m_runId == 0) {
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if (mMultiModal::getPopInfor())
			mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
#endif // OFEC_CONSOLE

	cmaes_t evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xfinal;
	int i;
	
	int numDim = Global::msp_global->mp_problem->getNumDim();

	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(&evo, numDim, NULL, NULL, 0, m_pop.size(), initialsFilePathName.c_str());//"cmaes_initials.par"
	//printf("%s\n", cmaes_SayHello(&evo));
	cmaes_ReadSignals(&evo, signalsFilePathName.c_str());  /* "cmaes_signals.par"write header and initial values */

												   /* Iterate until stop criterion holds */
	while (!ifTerminating())/*!cmaes_TestForTermination(&evo)*/
	{
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */
		for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
			copy(pop[i], m_pop[i].data());
		}

		for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) {
			/* You may resample the solution i until it lies within the
			feasible domain here, e.g. until it satisfies given
			box constraints (variable boundaries). The function
			is_feasible() needs to be user-defined.
			Assumptions: the feasible domain is convex, the optimum
			is not on (or very close to) the domain boundary,
			initialX is feasible (or in case typicalX +- 2*initialStandardDeviations
			is feasible) and initialStandardDeviations is (are)
			sufficiently small to prevent quasi-infinite looping.
			*/
			while (!Global::msp_global->mp_problem->isValid(m_pop[i].data())) {
				cmaes_ReSampleSingle(&evo, i);
				copy(pop[i], m_pop[i].data());
			}		
		}
											
		/* evaluate the new search points using fitfun */
		for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
			arFunvals[i] = fitCompute(m_pop[i]);
		}

		/* update the search distribution used for cmaes_SamplePopulation() */
		cmaes_UpdateDistribution(&evo, arFunvals);

		/* read instructions for printing output or changing termination conditions */
		cmaes_ReadSignals(&evo, signalsFilePathName.c_str());//"cmaes_signals.par"
		fflush(stdout); /* useful in MinGW */

#ifdef OFEC_CONSOLE
		if (mMultiModal::getPopInfor()) {
			int peaksf;
			peaksf = CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
			mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(), \
				Global::msp_global->m_totalNumIndis, 1, peaksf, \
				CAST_PROBLEM_CONT->getGOpt().getNumOpt(), 0, 0, 0, 0, \
				0, 0, CAST_PROBLEM_CONT->getGOpt().isAllFound());
		}
#endif
		//cout << cmaes_Get(&evo, "fbestever") <<" "<<Global::msp_global->mp_problem->getEvaluations()<< endl;
	}
	//printf("Stop:\n%s\n", cmaes_TestForTermination(&evo)); /* print termination reason */
	//cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

															/* get best estimator for the optimum, xmean */
	xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
	cmaes_exit(&evo); /* release memory */

					  /* do something with final solution and finally release memory */
	free(xfinal);

    return Return_Terminate;
}
