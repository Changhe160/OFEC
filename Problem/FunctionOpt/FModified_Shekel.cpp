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
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FModified_Shekel.h"

FModified_Shekel::FModified_Shekel(ParamMap &v):Problem(v[param_proId], v[param_numDim],v[param_proName],1),\
	BenchmarkFunction(v[param_proId], v[param_numDim],v[param_proName],1){
	if(m_numDim>5) throw myException("number of dim must be <=5@ FModified_Shekel::FModified_Shekel");
	 setSearchRange(0.0,11.0);
     initialize();
}
FModified_Shekel::FModified_Shekel(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
	if(m_numDim>5) throw myException("number of dim must be <=5@ FModified_Shekel::FModified_Shekel");
	 setSearchRange(0.0,11.0);
     initialize();
}

FModified_Shekel::~FModified_Shekel(){
    //dtor
}
void FModified_Shekel::initialize(){
	double a[8][5]={4,4,6.3,4,4,1,1,8.5,1,1,6,6,9.1,6,6,3.5,7.5,4,9,4,5,5,3,3,9,9.1,8.2,2,3,9,1.5,9.3,7.4,3,9,7.8,2.2,5.3,9,3};
	double c[8]={0.1,0.2,0.4,0.15,0.6,0.2,0.06,0.18};
	copy(c,c+8,m_c);
	for(int i=0;i<8;i++) copy(a[i],a[i]+5,m_a[i]);
	setOptType(MAX_OPT);
	setDisAccuracy(0.2);
	setAccuracy(1.e-3);
	m_globalOpt.setFlagLocTrue();	
	m_globalOpt.setNumOpts(8); //1 gopt+7 lopt
	 
	CodeVReal x(m_numDim,1);
	for(int i=0;i<m_numDim;++i) x.m_x[i]=m_a[6][i];
	BenchmarkFunction::evaluate_(x,false);
	m_globalOpt[0].data()=x;
	for(int m=0,k=0;m<8;m++){
		if(m==6) continue;
		for(int i=0;i<m_numDim;++i) x.m_x[i]=m_a[m][i];
		BenchmarkFunction::evaluate_(x,false);
		m_globalOpt[++k].data()=x;
	}
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FModified_Shekel::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	for(int i=0;i<8;++i){
		double b=0;
		for(int j=0;j<m_numDim;++j){
			b+=(x[j]-m_a[i][j])*(x[j]-m_a[i][j]);
		}
		s+=pow(b+m_c[i],-1);
	}
	
	obj[0]= s+m_bias;
}