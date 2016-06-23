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

#ifndef NRAND1DEINDIVIDUAL_H
#define NRAND1DEINDIVIDUAL_H

#include "../DEIndividual.h"

class Nrand1DEIndividual :public DEIndividual
{
public:
	Nrand1DEIndividual();
	~Nrand1DEIndividual(){}
	ReturnFlag run_();
	void recombine(double CR);
};

#endif //NRAND1DEINDIVIDUAL_H