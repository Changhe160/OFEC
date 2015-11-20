#include "NSGAII_SBXRealMu.h"


NSGAII_SBXRealMu::NSGAII_SBXRealMu(ParamMap &v): NSGAII<GAIndividual<CodeVReal>,SBX_RealMuPopulation>(v)
{
	m_parent.setCrossXP(0.9);
	m_parent.setEta(20,20);
}

void NSGAII_SBXRealMu::evolve_mo()
{
	int m=0;
	if(m_parent.getPopSize()%2!=0)
		throw myException("population size should be even @NSGAII_SBXRealMu::evolve_mo()");
	for(int n=0; n<m_parent.getPopSize(); n+=2)
	{
		vector<int> p(2);       
		p[0] = tour_selection();
		while(1){  p[1] = tour_selection();  	if(p[1]!=p[0]) break; }

		m_parent.cross_mutate(p,m_offspring[m],m_offspring[m+1]);
		m+=2;
		*m_offspring[m++]=*m_parent[n];
		*m_offspring[m++]=*m_parent[n+1];
	}
}

