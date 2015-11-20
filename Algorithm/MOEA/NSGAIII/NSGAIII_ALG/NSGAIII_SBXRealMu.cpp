#include "NSGAIII_SBXRealMu.h"


NSGAIII_SBXRealMu::NSGAIII_SBXRealMu(ParamMap &v): NSGAIII<GAIndividual<CodeVReal>,SBX_RealMuPopulation>(v)
{

}

void NSGAIII_SBXRealMu::evolve_mo()
{
	int m=m_parent.getPopSize();
	if(m_parent.getPopSize()%2!=0)
		throw myException("population size should be even @NSGAIII_SBXRealMu::evolve_mo()");
	for (size_t n=0; n<m_parent.getPopSize(); n+=2)
	{
		vector<int> p(2); 
		p[0] = Global::msp_global->getRandInt(0,m_parent.getPopSize());
		p[1] = Global::msp_global->getRandInt(0,m_parent.getPopSize());

		m_parent.cross_mutate(p,m_offspring[m],m_offspring[m+1]);
		m+=2;
	}
	for(int i=0;i<m_parent.getPopSize();i++)
		*m_offspring[i]=*m_parent[i];
}

