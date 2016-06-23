
#ifndef FAMO_SWARM__H
#define FAMO_SWARM__H

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