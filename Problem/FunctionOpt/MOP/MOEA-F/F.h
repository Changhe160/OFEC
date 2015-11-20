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

/*==========================================================================
//  Implementation of MOEA/D Based on Differential Evolution (DE) for Continuous Multiobjective 
//  Optimization Problems with Complicate Pareto Sets (2007)
//
//  See the details of MOEA/D-DE and test problems in the following paper
//  H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization 
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  The component functions of each test instance can be found in "objective.h". 
//
//  The source code of MOEA/D-DE and NSGA-II-DE were implemented by Hui Li and Qingfu Zhang  
//
//  If you have any questions about the codes, please contact 
//  Qingfu Zhang at qzhang@essex.ac.uk  or Hui Li at hzl@cs.nott.ac.uk
===========================================================================*/

#ifndef FBASE_H
#define FBASE_H

#include "../../BenchmarkFunction.h"


class F_Base :public BenchmarkFunction
{
public:
	F_Base(int ID, int numDim, const string &proName, int numObj);
	~F_Base(){}
	void evaluate__(double const *x,vector<double>& obj);
	int getDtype() const { return m_dtype; }
	int getPtype() const { return m_ptype; }
	int getLtype() const { return m_ltype; }
protected:
	void alphafunction(double alpha[], double const *x, int dim, int type);
	double betafunction(const vector<double> & x, int type);
	double psfunc2(const double &x,const double &t1, int dim, int type, int css);
	double psfunc3(const double &x,const double &t1,const double &t2, int dim, int type);
	void calObjective(double const *x_var, vector <double> &y_obj);
	void LoadPF();

	int m_dtype, m_ptype, m_ltype;
};

#endif //FBASE_H