#ifndef CPSOSUBSWARM_H
#define CPSOSUBSWARM_H
#include "../../Particle.h"
#include "../../Swarm.h"

class CPSOSwarm;
class CPSOSubSwarm;
class CPSOParticle:public Particle{
	friend class CPSOSubSwarm;
	friend class CPSOSwarm;
public:
	CPSOParticle();
	~CPSOParticle();
	CPSOParticle( const Solution<CodeVReal> &chr);
};

class CPSOSubSwarm:public Swarm<CodeVReal,CPSOParticle> {	
public:
	friend class CPSOSwarm;
	CPSOSubSwarm(int size, bool flag);
	~CPSOSubSwarm(){}
	CPSOSubSwarm(Group<CodeVReal,CPSOParticle> &g);
protected:
	ReturnFlag evolve(); 
};

#endif