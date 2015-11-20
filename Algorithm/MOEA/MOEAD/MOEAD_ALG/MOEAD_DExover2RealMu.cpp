#include "MOEAD_DExover2RealMu.h"

MOEAD_DExover2RealMu::MOEAD_DExover2RealMu(ParamMap &v) :MOEAD<DEIndividual,DExover2_RealMuPopulation>(v)
{
	m_best.clear();
	m_bestIdx.clear();
	m_worstIdx.clear();
}

void MOEAD_DExover2RealMu::evolve_mo()
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
		matingselection(p,n,2,type);  // neighborhood selection

		// produce a child TypeIndiv
		DEIndividual child;
		vector<int> index(3);
		index[0]=n; index[1]=p[0]; index[2]=p[1];
		this->cross_mutate(index,child);

		// update the reference points and other TypeIndivs in the neighborhood or the whole population
		update_reference(child);
		update_problem(child, n, type);
	}
} 

