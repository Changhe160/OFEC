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
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FIBA.h"


FIBA::FIBA(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1),m_case(v[param_case]=(v.find(param_case)==v.end()?1:(int)v[param_case])){
	 v[param_numDim]=2;
	 setSearchRange(-4.0,4.0);
     initialize();
}
FIBA::FIBA(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1),m_case(1){
	 setSearchRange(-4.0,4.0);
     initialize();
}
void FIBA::setPara(){
	if(m_case==1){
		m_kappa=-0.95;
		m_chi=-1.26;
	}else{
		m_kappa=4;
		m_chi=2;
	}
}
void FIBA::setCase(int c){
	m_case=c;
	setPara();
}
FIBA::~FIBA(){
    //dtor
}
void FIBA::initialize(){
	setPara();
	m_globalOpt.setFlagLocTrue();	
	m_globalOpt.setNumOpts(4);
	setAccuracy(1.e-6);
	CodeVReal x(m_numDim,1);
	if(m_case==1){	
		setDisAccuracy(0.08);
		x.m_x[0]=0.45186;x.m_x[1]=0.0; x.m_obj[0]=-0.0080161;
		m_globalOpt[0].data()=x;
		x.m_x[0]=-0.22593;x.m_x[1]=0.39132; 
		m_globalOpt[1].data()=x;
		x.m_x[0]=-0.22593;x.m_x[1]=-0.39132; 
		m_globalOpt[2].data()=x;
		x.m_x[0]=0.0;x.m_x[1]=0.0; x.m_obj[0]=0;
		m_globalOpt[3].data()=x;
	}else{
		setDisAccuracy(0.5);
		x.m_x[0]=0.0;x.m_x[1]=0.0; x.m_obj[0]=0;
		m_globalOpt[0].data()=x;
		x.m_x[0]=1.2243;x.m_x[1]=0.0; x.m_obj[0]=0.71449;
		m_globalOpt[1].data()=x;
		x.m_x[0]=-0.61215;x.m_x[1]=1.0603; 
		m_globalOpt[2].data()=x;
		x.m_x[0]=-0.61215;x.m_x[1]=-1.0603; 
		m_globalOpt[3].data()=x;
	}
	addProTag(MMP);
	setObjSet();
}
void FIBA::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	double t0=x[0]*x[0]+x[1]*x[1];
	s=(t0)/(1+t0)+m_kappa*(14*t0+pow(t0,2)*m_chi*m_chi-2*sqrt(14.)*(pow(x[0],3)-3*x[0]*x[1]*x[1])*m_chi)/(14*(pow(1+t0,2)));
	obj[0]= s+m_bias;

}