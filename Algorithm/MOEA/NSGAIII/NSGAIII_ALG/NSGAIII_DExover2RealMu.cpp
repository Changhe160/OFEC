#include "NSGAIII_DExover2RealMu.h"

NSGAIII_DExover2RealMu::NSGAIII_DExover2RealMu(ParamMap &v) :NSGAIII<DEIndividual,DExover2_RealMuPopulation>(v)
{
}

void NSGAIII_DExover2RealMu::evolve_mo()
{
	int m=0;
	for(int n=0; n<m_parent.getPopSize(); n++)
	{
		vector<int> p(3);       
		p[0] = Global::msp_global->getRandInt(0,m_parent.getPopSize());
		p[1] = Global::msp_global->getRandInt(0,m_parent.getPopSize());
		p[2] = Global::msp_global->getRandInt(0,m_parent.getPopSize());

		m_parent.cross_mutate(p,*m_offspring[m]);
		++m;
		*m_offspring[m++]=*m_parent[n];
	}
} 