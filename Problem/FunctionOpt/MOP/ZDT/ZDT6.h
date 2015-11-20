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

#ifndef ZDT6_H
#define ZDT6_H


#include "ZDT.h"

class ZDT6 :public ZDT
{
public:
	ZDT6(ParamMap &v);
	~ZDT6(){};
protected:
	void evaluate__(double const *x,vector<double>& obj);
};

#endif //ZDT6_H