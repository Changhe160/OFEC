/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 13 Jan 2015
// Last modified:

#ifndef MOEAD_SBX_REALMU_H
#define MOEAD_SBX_REALMU_H

#include "../MOEAD.h"
#include "../../../GA/MOEA_GA/SBX_RealMuPopulation.h"

class MOEAD_SBXRealMu :public MOEAD<GAIndividual<CodeVReal>,SBX_RealMuPopulation>
{
public:
	MOEAD_SBXRealMu(ParamMap &v);
	~MOEAD_SBXRealMu(){};
protected:
	void evolve_mo();	
};


#endif //MOEAD_SBX_REALMU_H