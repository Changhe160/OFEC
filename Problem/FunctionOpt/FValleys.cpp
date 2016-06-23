/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FValleys.h"


FValleys::FValleys(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	 vector<double> lower, upper;
	  v[param_numDim]=2;
     upper.push_back (3);
	 upper.push_back (2);
	 lower.push_back (-2.5);
	 lower.push_back (-2);
	 setSearchRange(lower,upper);
     initialize();
}
FValleys::FValleys(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	vector<double> lower, upper;
	 upper.push_back (3);
	 upper.push_back (2);
	 lower.push_back (-2.5);
	 lower.push_back (-2);
	 setSearchRange(lower,upper);
     initialize();
}

FValleys::~FValleys(){
    //dtor
}
void FValleys::initialize(){
	setOptType(MAX_OPT); 
	setAccuracy(1.e-4);
	setDisAccuracy(0.5);
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(2); //1 gopt + 1 lopt
	 
	CodeVReal x(m_numDim,1);
	x.m_x[0]=1.69714;x.m_x[1]=0.0; x.m_obj[0]=4.8168;
	m_globalOpt[0].data()=x;
	x.m_x[0]=-1.44446;x.m_x[1]=0.0; x.m_obj[0]=3.2460;
	m_globalOpt[1].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FValleys::evaluate__(double const *x,vector<double>& obj){
	double s;
	s=sin(2*x[0]-0.5*OFEC_PI)+3*cos(x[1])+0.5*x[0];
	obj[0]= s;

}