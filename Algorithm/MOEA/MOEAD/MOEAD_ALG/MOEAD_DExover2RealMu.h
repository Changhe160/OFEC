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

#ifndef MOEAD_DEXOVER2REALMU_H
#define MOEAD_DEXOVER2REALMU_H

#include "../MOEAD.h"
#include "../../../DE/MOEA_DE/DExover2_RealMuPopulation.h"

class MOEAD_DExover2RealMu :public MOEAD<DEIndividual,DExover2_RealMuPopulation>
{
public:
	MOEAD_DExover2RealMu(ParamMap &v);
	~MOEAD_DExover2RealMu(){};
protected:
	void evolve_mo();
};


#endif //MOEAD_DEXOVER2REALMU_H
