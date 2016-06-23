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
// Created: 31 December 2014
// Last modified:

#ifndef ZDT1_H
#define ZDT1_H


#include "ZDT.h"

class ZDT1 :public ZDT
{
public:
	ZDT1(ParamMap &v);
	~ZDT1(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //ZDT1_H