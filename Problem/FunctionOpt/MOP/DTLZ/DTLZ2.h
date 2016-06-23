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

#ifndef DTLZ2_H
#define DTLZ2_H

#include "DTLZ.h"

class DTLZ2 :public DTLZ
{
public:
	DTLZ2(ParamMap &v);
	~DTLZ2(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //DTLZ2_H                                        