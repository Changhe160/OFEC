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

#ifndef NSGAII_SBX_REALMU_H
#define NSGAII_SBX_REALMU_H

#include "../NSGAII.h"
#include "../../../GA/MOEA_GA/SBX_RealMuPopulation.h"

class NSGAII_SBXRealMu :public NSGAII<GAIndividual<CodeVReal>,SBX_RealMuPopulation>
{
public:
	NSGAII_SBXRealMu(ParamMap &v);
	~NSGAII_SBXRealMu(){};
protected:
	void evolve_mo();	
};


#endif //NSGAII_SBX_REALMU_H