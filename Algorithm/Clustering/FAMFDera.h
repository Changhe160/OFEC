#ifndef FAMF_DERATING
#define FAMF_DERATING

#include "../../Problem/ContOptima.h"
#include "../../Problem/optimum.h"
struct FAMFDerating
{
	#ifdef OFEC_CONSOLE
	static boost::thread_specific_ptr<Optima<CodeVReal,ContOptimum>> msp_opt;
	#endif
	#ifdef OFEC_DEMON
	static unique_ptr<Optima<CodeVReal,ContOptimum>> msp_opt;
	#endif
	
	static bool ms_enableDerating;
	static double getDeratingFactor(double dis, double ref);
	template<typename ED>
	static bool derateFitness(Solution<ED> & x);
};

inline double FAMFDerating::getDeratingFactor(double dis, double ref){
	if(ref<=0) return 0;
	double f=dis*dis*dis/(ref*ref*ref);
	f=exp(-f)*(1-f);

	return f>0?f:0;;
}
template<typename ED>
bool FAMFDerating::derateFitness(Solution<ED> & x){

	if(!ms_enableDerating) return false;
	if(msp_opt->getNumOpt()<=0) return false;


	int idx;
	double dis;
	msp_opt->getNearestGOpt(x,&idx,&dis);
	double refdis=(*msp_opt)[idx].getRefRadi(x);
	if(dis<refdis){
		x.data().m_obj=(*msp_opt)[idx].getDeratingObj();
		return true;
	}
	return false;
}

#endif