/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 21 September 2011
// Last modified: 21 Nov 2014

/*
SPSO07: C. Maurice, "Standard pso 2007 (spso-07)" http://www.particleswarm.info/Programs.html, 2007.
*/
#ifndef SWARMLBEST_H
#define SWARMLBEST_H

#include "../Particle.h"
#include "../Swarm.h"

class SPSO07;
class ParticleLbest: public Particle{
	friend class SPSO07;
public:
	ReturnFlag move(double w, double c1, double c2, const Solution<CodeVReal> *gbest=0,bool clamping=false);
};

class SPSO07: public Swarm<CodeVReal,ParticleLbest>{
public:
    SPSO07();
    ~SPSO07();
	SPSO07(ParamMap&);
	ReturnFlag run_();           
protected:
	ReturnFlag evolve();
	void setNeibourhood();
	Solution<CodeVReal>& neighborBest(int idx);
private:
    unsigned m_impr;
};

#endif // SWARMLBEST_H
