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

/*
K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, "Scalable test problems
for evolutionary multi - objective optimization, " in Evolutionary Multiobjective
Optimization, A.Abraham, L.Jain, and R.Goldberg, Eds.
London, U.K.: Springer - Verlag, 2005, pp. 105¨C145.
*/

#ifndef DTLZ_H
#define DTLZ_H

#include "../../BenchmarkFunction.h"

class DTLZ :public BenchmarkFunction
{
public:
	DTLZ(int ID, int numDim, const string &proName, int numObj);
	~DTLZ(){};
protected:
	void evaluate__(double const *x,vector<double>& obj){}
	void generateAdLoadPF();
};

#endif //DTLZ_H