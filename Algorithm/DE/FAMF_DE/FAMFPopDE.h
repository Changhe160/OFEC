#ifndef FAMF_POP_DE_H
#define FAMF_POP_DE_H

#include "FAMFIndiDE.h"
#include "../DEPopulation.h"
#include "../../Clustering/FAMFPop.h"

template<typename, typename, typename> class FAMF;

class FAMFPopDE: public FAMFPop<FAMFIndiDE,DEPopulation<CodeVReal,FAMFIndiDE>>{	
	template<typename, typename, typename> friend class FAMF;
public:	
	FAMFPopDE(Group<CodeVReal,FAMFIndiDE> &g):FAMFPop<FAMFIndiDE,DEPopulation<CodeVReal,FAMFIndiDE>>(g){
		setMutationStrategy(DE_best_2); 
		setParmeter(0.6,0.5);
	}
protected:
	ReturnFlag evolve();
};
#endif