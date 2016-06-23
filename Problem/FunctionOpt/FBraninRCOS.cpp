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
// Created: 26 Dec 2014
// Last modified:

#include "FBraninRCOS.h"


FBraninRCOS::FBraninRCOS(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	v[param_numDim]=2;
	vector<double> lower, upper;
	upper.push_back (10);
	upper.push_back (15);
	lower.push_back (-5);
	lower.push_back (0);
	setSearchRange(lower,upper);
	initialize();
}
FBraninRCOS::FBraninRCOS(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	vector<double> lower, upper;
	upper.push_back (10);
	upper.push_back (15);
	lower.push_back (-5);
	lower.push_back (0);
	setSearchRange(lower,upper);
	initialize();
}

FBraninRCOS::~FBraninRCOS(){
    //dtor
}
void FBraninRCOS::initialize(){
	setOptType(MIN_OPT);
	setDisAccuracy(1.0);
	setAccuracy(1.e-5);
	m_globalOpt.setFlagLocTrue();	 
	m_globalOpt.setNumOpts(3);
	
	CodeVReal x(m_numDim,1);
	x.m_x[0]=-OFEC_PI;x.m_x[1]=12.275; x.m_obj[0]=0.397887;
	evaluate_(x,false);
	m_globalOpt[0].data()=x;
	x.m_x[0]=OFEC_PI;x.m_x[1]=2.275; 
	m_globalOpt[1].data()=x;
	x.m_x[0]=9.42478;x.m_x[1]=2.475; 
	m_globalOpt[2].data()=x;
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FBraninRCOS::evaluate__(double const *x,vector<double>& obj){
	double s,a;
	a=x[1]-5.1*x[0]*x[0]/(4*OFEC_PI*OFEC_PI)+5*x[0]/OFEC_PI-6;
	s=a*a+10*(1-1/(8*OFEC_PI))*cos(x[0])+10;
	obj[0]= s+m_bias;

}