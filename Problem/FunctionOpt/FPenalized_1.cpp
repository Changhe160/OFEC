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
*************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FPenalized_1.h"

FPenalized_1::FPenalized_1(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-50,50);

     initialize();
}
FPenalized_1::FPenalized_1(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-50,50);

     initialize();
}

FPenalized_1::~FPenalized_1(){
    //dtor
}


void FPenalized_1::initialize(){
	vector<double> v(m_numDim,-1);
    setOriginalGlobalOpt(0,&v);
    setGlobalOpt(0,&v);
}

void FPenalized_1::evaluate__(double const *x,vector<double>& obj){
    double * y = new double[m_numDim];
	for(int i=0;i<m_numDim;i++) y[i]=(x[i]+1)/4.+1;
	double s=0;
	for(int i=0;i<m_numDim-1;i++)
		s+=(y[i]-1)*(y[i]-1)*(1+10*sin(OFEC_PI*y[i+1])*sin(OFEC_PI*y[i+1]));
	s+=(y[m_numDim-1]-1)*(y[m_numDim-1]-1)+10*sin(OFEC_PI*y[0])*sin(OFEC_PI*y[0]);
	s=s*OFEC_PI/m_numDim;
	for(int i=0;i<m_numDim;i++) {
		s+=u(x[i],10,100,4);
	}
	delete []y;
	y=0;
	obj[0]= s+m_bias;

}
double FPenalized_1::u(double x, double a, double k, double m)const{
	if(x>a) return k*pow(x-a,m);
	else if(x<-a) return k*pow(-x-a,m);
	else return 0;
}
