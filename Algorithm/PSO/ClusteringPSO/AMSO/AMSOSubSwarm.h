#ifndef AMSOSUBSWARM_H
#define AMSOSUBSWARM_H
#include "../../Particle.h"
#include "../../Swarm.h"

class AMSOSwarm;
class AMSOSubSwarm;
class AMSOParticle:public Particle{
	friend class AMSOSubSwarm;
	friend class AMSOSwarm;
public:
	AMSOParticle();
	AMSOParticle( const Solution<CodeVReal> &chr);
};

class AMSOSubSwarm:public Swarm<CodeVReal,AMSOParticle> {	
public:
	friend class AMSOSwarm;
	AMSOSubSwarm(int size, bool flag);
	~AMSOSubSwarm(){}
	AMSOSubSwarm(Group<CodeVReal,AMSOParticle> &g);
protected:
	ReturnFlag evolve(); 
};

#endif