#include "NSGAII_DExover2RealMu.h"

NSGAII_DExover2RealMu::NSGAII_DExover2RealMu(ParamMap &v) :NSGAII<DEIndividual,DExover2_RealMuPopulation>(v)
{

}

void NSGAII_DExover2RealMu::evolve_mo()
{
	int m=0;
	for(int n=0; n<m_parent.getPopSize(); n++)
	{
		vector<int> p(3);       
		p[0] = tour_selection();
		while(1){  p[1] = tour_selection();  	if(p[1]!=p[0]) break; }
		while(1){  p[2] = tour_selection();  	if(p[2]!=p[0]&&p[2]!=p[1]) break; }

		m_parent.cross_mutate(p,*m_offspring[m]);
		++m;
		*m_offspring[m++]=*m_parent[n];
	}
} 
