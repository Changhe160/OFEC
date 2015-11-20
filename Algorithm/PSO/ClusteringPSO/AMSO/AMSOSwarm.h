#ifndef AMSOSWARM_H
#define AMSOSWARM_H

#include "AMSOSubSwarm.h"
#include "../../../../Utility/Clustering/HSLHClustering.h"
#include "../../../MultiPopulationCont.h"

/*
C. Li, S. Yang and M. Yang, "An Adaptive Multi-Swarm Optimizer for Dynamic
Optimization Problems," MIT.
*/

class AMSOSwarm: public Swarm<CodeVReal,AMSOParticle>, public MultiPopulationCont<AMSOSubSwarm> {
	struct PopInfor{
	// record status of the whole populations, including current number of fitness evalutions, number of pops
	int m_evals;
	int m_pops;
	PopInfor(int evals,int pops):m_evals(evals),m_pops(pops){}
	PopInfor(const PopInfor & infor){m_evals=infor.m_evals;m_pops=infor.m_pops;}
	PopInfor & operator=(const PopInfor & infor){
		if(this==&infor) return *this;
		m_evals=infor.m_evals;m_pops=infor.m_pops;
		return *this;
	}
	};
private:
	queue<PopInfor>mq_pops;	// information of whole populations during one diversity increase
	int m_maxIndis,m_minIndis,m_stepIndis,m_gap;
	//maximum, minimum number of indis allowed in the search space;
	//m_stepIndis: number of indis to be added or removed when necessary;
	//m_gap: number of fitness evalutations between two successive diversity increase operation
    const double m_alpha;  // mc_alpha: threshold value for the convergence of pop numbers;
	const int m_offPeak; //threshold of the difference of pop numbers between two successive diversity increase

	int m_convergedPops,m_preIndis,m_nextIndis;
	// total individual at prevous and next diversity increase

	int m_bufferTimes;  //counter of generations where there is no population change in two successive check points.
	int m_maxPops;    //largest number of pops since last diversity increase
	vector<unique_ptr<AMSOParticle>> mv_converge;

	Cluster<CodeVReal,AMSOParticle> m_clst;
	int m_initialSize,m_subSize;
public:
	AMSOSwarm(ParamMap &v);
	~AMSOSwarm();
	ReturnFlag run_();
private:
	void initialize();
	void createSubswarms();
	void getNextIndis(PopInfor &infor);
	void increaseDiversity();
};

#endif