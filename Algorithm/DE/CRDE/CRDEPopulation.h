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
// Created: 15 Jan 2015
// Last modified:
/*
Thomsen R (2004) Multimodal optimization using crowding-based differential evolution. In: Proceedings of
the 2004 IEEE congress on evolutionary computation, pp 1382¨C1389. IEEE Press, Portland, OR, USA
*/
#ifndef CRDEPOPULATION_H
#define CRDEPOPULATION_H

#include "../DEIndividual.h"
#include "../DEPopulation.h"

class CRDEPopulation :public DEPopulation<CodeVReal,DEIndividual>
{
protected:
	void updateDis(int i);
public:
	CRDEPopulation(ParamMap &v);
	~CRDEPopulation(){}
    ReturnFlag run_();
    ReturnFlag evolve();
	bool ifTerminating();
};

#endif // CRDEPOPULATION_H
