#include "MOEAD_SBXRealMu.h"


MOEAD_SBXRealMu::MOEAD_SBXRealMu(ParamMap &v): MOEAD<GAIndividual<CodeVReal>,SBX_RealMuPopulation>(v)
{
	m_best.clear();
	m_bestIdx.clear();
	m_worstIdx.clear();
}

void MOEAD_SBXRealMu::evolve_mo()
{
	vector<int> perm(this->m_popsize);
	Global::msp_global->initializeRandomArray<vector<int> >(perm,this->m_popsize);

	for(int i=0; i<this->m_popsize; i++)
	{
		int n = perm[i];
		// or int n = i;
		int type;
		double rnd = Global::msp_global->mp_uniformAlg->Next();

		// mating selection based on probability
		if(rnd<m_realb)    type = 1;   // neighborhood
		else             type = 2;   // whole population

		// select the indexes of mating parents
		vector<int> p;
		matingselection(p,n,1,type);  // neighborhood selection

		// produce a child TypeIndiv
		GAIndividual<CodeVReal> child1,child2;
		vector<int> index(2);
		index[0]=n; index[1]=p[0];
		this->cross_mutate(index,&child1,&child2);

		// update the reference points and other TypeIndivs in the neighborhood or the whole population
		update_reference(child1);
		update_reference(child2);
		update_problem(child1, n, type);
		update_problem(child2, n, type);
	}
}

