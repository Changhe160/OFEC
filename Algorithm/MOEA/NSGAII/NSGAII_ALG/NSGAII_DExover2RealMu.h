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

#ifndef NSGAII_DEXOVER2REALMU_H
#define NSGAII_DEXOVER2REALMU_H

#include "../NSGAII.h"
#include "../../../DE/MOEA_DE/DExover2_RealMuPopulation.h"

class NSGAII_DExover2RealMu :public NSGAII<DEIndividual,DExover2_RealMuPopulation>
{
public:
	NSGAII_DExover2RealMu(ParamMap &v);
	~NSGAII_DExover2RealMu(){};
protected:
	void evolve_mo();	
};


#endif //NSGAII_DEXOVER2REALMU_H
