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

#ifndef DTLZ4_H
#define DTLZ4_H

#include "DTLZ.h"

class DTLZ4 :public DTLZ
{
public:
	DTLZ4(ParamMap &v);
	~DTLZ4(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //DTLZ4_H                                        