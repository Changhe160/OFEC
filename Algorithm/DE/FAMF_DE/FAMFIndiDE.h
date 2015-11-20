#ifndef FAMF_INDI_DE_H
#define FAMF_INDI_DE_H
#include "../DEIndividual.h"

class FAMFPopDE;
template<typename, typename, typename> class FAMF;
class FAMFIndiDE: public DEIndividual{
friend class FAMFPopDE;
template<typename, typename, typename> friend class FAMF;
public:
	FAMFIndiDE();
	FAMFIndiDE( const Solution<CodeVReal> &chr);
	// Brownian movements to dealwith noisy environments for the gbest particle
	ReturnFlag brwonianMove(double radius);
	ReturnFlag select();
	ReturnFlag cauchyMove(double radius=-1);
};
#endif