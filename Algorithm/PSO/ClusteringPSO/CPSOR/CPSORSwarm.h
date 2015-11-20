#ifndef CPSORSWARM_H
#define CPSORSWARM_H

#include "../../../../Utility/Clustering/HSLHClustering.h"
#include "../../../MultiPopulationCont.h"
#include "CPSORSubSwarm.h"

/*
C. Li and S. Yang, "A general framework of multi-population methods with clustering in undetectable dynamic environments,"
IEEE Trans. Evol. Comput., vol. 16, no. 4, pp. 556¨C577, 2012.
*/

class CPSORSwarm: public Swarm<CodeVReal,CPSORParticle>, public MultiPopulationCont<CPSORSubSwarm> {
private:
	Cluster<CodeVReal,CPSORParticle> m_clst;
	int m_initialSize,m_subSize;
	double m_diversityDegree;
public:
	CPSORSwarm(ParamMap &v);
	~CPSORSwarm(){};
	ReturnFlag run_();
private:
	void initialize();
	void createSubswarms( );
};

#endif