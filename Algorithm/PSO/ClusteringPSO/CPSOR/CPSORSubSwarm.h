#ifndef CPSORSUBSWARM_H
#define CPSORSUBSWARM_H
#include "../../Particle.h"
#include "../../Swarm.h"

class CPSORSwarm;
class CPSORSubSwarm;
class CPSORParticle:public Particle{
	friend class CPSORSubSwarm;
	friend class CPSORSwarm;
public:
	CPSORParticle();
	CPSORParticle( const Solution<CodeVReal> &chr);
};

class CPSORSubSwarm:public Swarm<CodeVReal,CPSORParticle> {	
public:
	friend class CPSORSwarm;
	CPSORSubSwarm(int size, bool flag);
	~CPSORSubSwarm(){}
	CPSORSubSwarm(Group<CodeVReal,CPSORParticle> &g);
protected:
	ReturnFlag evolve(); 
};

#endif