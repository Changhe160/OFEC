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
// Created: 14 Jan 2015
// Last modified:

#ifndef NRAND1DEPOPULATION_H
#define NRAND1DEPOPULATION_H
/*
M. Epitropakis, V. Plagianakos, and M. Vrahatis, ¡°Finding multiple
global optima exploiting differential evolution¡¯s niching capability,¡± in
IEEE Symposium on Differential Evolution, 2011. SDE 2011. (IEEE
Symposium Series on Computational Intelligence), Paris, France, April
2011, p. 80¨C87.
*/
#include "../DEPopulation.h"
#include "Nrand1DEIndividual.h"

class Nrand1DEPopulation :public DEPopulation<CodeVReal,Nrand1DEIndividual>
{
public:
	Nrand1DEPopulation(ParamMap &v);
	~Nrand1DEPopulation(){}
	ReturnFlag run_();
	bool ifTerminating();
protected:
	void mutate(const int idx);

};


#endif //NRAND1DEPOPULATION_H