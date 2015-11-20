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
// Created: 31 December 2014
// Last modified:

#ifndef F3_H
#define F3_H

#include "F.h"


class F3 :public F_Base
{
public:
	F3(ParamMap &v);
	~F3(){}
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //F3_H