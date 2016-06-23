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
SPSO11: C. Maurice, "Standard PSO 2011 (SPSO-2011)" http://www.particleswarm.info/Programs.html, 2011.
*/
#ifndef SPSO11_H
#define SPSO11_H

#include "../Particle.h"
#include "../Swarm.h"

//M. Clerc, Standard Particle Swarm Optimisation,¡± Particle Swarm Central, Tech.Rep., 2012,
class SPSO11;
class ParticleSPSO11 : public Particle{
	friend class SPSO11;
public:
	void initializeVelocity();
	ReturnFlag move(double w, double c1, double c2, const Solution<CodeVReal> *lbest=0,bool clamping=false);
};

class SPSO11: public Swarm<CodeVReal, ParticleSPSO11>{
public:
    SPSO11();
	SPSO11(ParamMap&);
	SPSO11(int size) :SPSO11(Global::g_arg) {}
	ReturnFlag run_();           
protected:
	ReturnFlag evolve();
	void setNeibourhood();
	Solution<CodeVReal>& neighborBest(int idx);
protected:
    unsigned m_impr;
	double m_p; //
};

#endif // SWARMLBEST_H
