
#ifndef FAMO_SWARM__H
#define FAMO_SWARM__H
/*
C. Li, T. T. Nguyen, M. Yang, M. Mavrovouniotis, and S. Yang. An adaptive multi-population framework for
locating and tracking multiple optima. IEEE Transactions on Evolutionary Computation, accepted 18 November 2015
*/
#include "FAMFParticle.h"
#include "../../Swarm.h"
#include "../../../Clustering/FAMFPop.h"

template<typename, typename, typename> class FAMF;

class FAMFPopPSO: public FAMFPop<FAMFParticle,Swarm<CodeVReal,FAMFParticle>>{	
	template<typename, typename, typename> friend class FAMF;
public:	
	FAMFPopPSO(Group<CodeVReal,FAMFParticle> &g);
	void initilizeVelAftClustering();
protected:
	ReturnFlag evolve();
};
#endif