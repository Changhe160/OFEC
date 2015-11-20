#ifndef CLUSTERINGPSO_H
#define CLUSTERINGPSO_H

#include "../../../../Utility/Clustering/HSLHClustering.h"
#include "../../../MultiPopulationCont.h"
#include "CPSOSubSwarm.h"

/*
S. Yang and C. Li, ¡°A clustering particle swarm optimizer for locating
and tracking multiple optima in dynamic environments,¡± IEEE Trans.
Evol. Comput., vol. 14, no. 6, pp. 959¨C974, 2010.
*/

class CPSOSwarm: public Algorithm, public MultiPopulationCont<CPSOSubSwarm>{
private:
	vector<unique_ptr<CPSOParticle>> mv_converge;
	Cluster<CodeVReal,CPSOParticle> m_clst;
	int m_initialSize,m_subSize;
public:
	CPSOSwarm(ParamMap &v);
	~CPSOSwarm();
	ReturnFlag run_();
	int removeOverlapping();
private:
	void initialize();
	void createSubswarms( );
};

#endif