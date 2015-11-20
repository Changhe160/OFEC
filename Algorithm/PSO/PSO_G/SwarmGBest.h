/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 16 Jan 2015
// Last modified:

#ifndef SWARMGBEST_H
#define SWARMGBEST_H

#include "../Particle.h"
#include "../Swarm.h"

class SwarmGBest :public Swarm<CodeVReal,Particle>
{
public:
    SwarmGBest();
	~SwarmGBest(){}
	SwarmGBest(ParamMap &v);
	SwarmGBest(int popsize, bool mode);
	SwarmGBest(const Solution<CodeVReal> & center, double radius, int rPopsize,bool mode);
	ReturnFlag run_(); 
	ReturnFlag evolve();
};

#endif //SWARMGBEST_H