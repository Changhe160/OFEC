#include "FAMFDera.h"

#ifdef OFEC_CONSOLE
boost::thread_specific_ptr<Optima<CodeVReal,ContOptimum>> FAMFDerating::msp_opt;
#endif
#ifdef OFEC_DEMON
unique_ptr<Optima<CodeVReal,ContOptimum>> FAMFDerating::msp_opt;
#endif

bool FAMFDerating::ms_enableDerating=false;


