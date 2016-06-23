/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 13 Jan 2015
// Last modified:

#ifndef NSGAIII_SBX_REALMU_H
#define NSGAIII_SBX_REALMU_H

#include "../NSGAIII.h"
#include "../../../GA/MOEA_GA/SBX_RealMuPopulation.h"

class NSGAIII_SBXRealMu :public NSGAIII<GAIndividual<CodeVReal>,SBX_RealMuPopulation>
{
public:
	NSGAIII_SBXRealMu(ParamMap &v);
	~NSGAIII_SBXRealMu(){};
protected:
	void evolve_mo();	
};


#endif //NSGAIII_SBX_REALMU_H