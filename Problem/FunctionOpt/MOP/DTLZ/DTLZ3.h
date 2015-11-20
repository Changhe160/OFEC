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
// Created: 12 JAN 2014
// Last modified:

#ifndef DTLZ3_H
#define DTLZ3_H

#include "DTLZ.h"

class DTLZ3 :public DTLZ
{
public:
	DTLZ3(ParamMap &v);
	~DTLZ3(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //DTLZ3_H                                        