/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
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
// Created: 15 May 2013
// Last modified:
#ifndef AMSOIIPARTICLE_H
#define AMSOIIPARTICLE_H

#include "../../Particle.h"

class FAMFPopPSO;
template<typename, typename, typename> class FAMF;
class FAMFParticle: public Particle{
friend class FAMFPopPSO;
template<typename, typename, typename> friend class FAMF;
public:
	FAMFParticle();
	FAMFParticle( const Solution<CodeVReal> &chr);
	void initializeVelocityAftClustering();
	// Brownian movements to dealwith noisy environments for the gbest particle
	ReturnFlag brwonianMove(double radius);
	ReturnFlag move( const Solution<CodeVReal> & lbest,  const Solution<CodeVReal> &gbest,double w, double c1, double c2);
	ReturnFlag cauchyMove(double radius=-1);
};

#endif