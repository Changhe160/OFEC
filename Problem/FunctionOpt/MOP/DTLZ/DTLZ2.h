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