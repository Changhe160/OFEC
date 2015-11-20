/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
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
#include "FSix_humpCamelBack.h"

FSix_humpCamelBack::FSix_humpCamelBack(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	v[param_numDim]=2;
	initialize();
}
FSix_humpCamelBack::FSix_humpCamelBack(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	
	initialize();
}

FSix_humpCamelBack::~FSix_humpCamelBack(){
    //dtor
}
void FSix_humpCamelBack::initialize(){
	vector<double> lower, upper;
	upper.push_back (1.9);
	upper.push_back (1.1);
	lower.push_back (-1.9);
	lower.push_back (-1.1);
	setSearchRange(lower,upper);

	setAccuracy(1.e-4);
	setDisAccuracy(0.1);
	
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(6);// 2gopt+ 4 lopt
	
	CodeVReal x(m_numDim,1);
	x.m_x[0]=-0.089842; x.m_x[1]=0.712656 ; x.m_obj[0]=-1.0316;
	m_globalOpt[0].data()=x;
	x.m_x[0]=0.089842; x.m_x[1]=-0.712656 ;  
	m_globalOpt[1].data()=x;
	x.m_x[0]=-1.70361;x.m_x[1]=0.796084; x.m_obj[0]=-0.21546;
	m_globalOpt[2].data()=x;
	x.m_x[0]=1.70361;x.m_x[1]=-0.796084; x.m_obj[0]=-0.21546;
	m_globalOpt[3].data()=x;
	x.m_x[0]=-1.6071;x.m_x[1]=-0.56865; x.m_obj[0]=2.10425;
	m_globalOpt[4].data()=x;
	x.m_x[0]=1.6071;x.m_x[1]=0.56865; x.m_obj[0]=2.10425;
	m_globalOpt[5].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
}
void FSix_humpCamelBack::evaluate__(double const *x,vector<double>& obj){
	double s=x[0]*x[0],t=x[1]*x[1];
	s=(4-2.1*s+pow(x[0],4)/3)*s+x[0]*x[1]+(-4+4*t)*t;
	obj[0]= s+m_bias;
}