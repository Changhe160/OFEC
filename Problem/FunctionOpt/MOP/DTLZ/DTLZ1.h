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
// Created: 12 JAN 2014
// Last modified:

#ifndef DTLZ1_H
#define DTLZ1_H

#include "DTLZ.h"

class DTLZ1 :public DTLZ
{
public:
	DTLZ1(ParamMap &v);
	~DTLZ1(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //DTLZ1_H                                        