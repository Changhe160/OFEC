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
*************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FPenalized_2.h"


FPenalized_2::FPenalized_2(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-50,50);

     initialize();
}
FPenalized_2::FPenalized_2(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-50,50);

     initialize();
}

FPenalized_2::~FPenalized_2(){
    //dtor
}


void FPenalized_2::initialize(){
	vector<double> v(m_numDim,1.0);
    setOriginalGlobalOpt(0,&v);
    setGlobalOpt(0,&v);
}


void FPenalized_2::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	for(int i=0;i<m_numDim-1;i++)
		s+=(x[i]-1)*(x[i]-1)*(1+sin(3*OFEC_PI*x[i+1])*sin(3*OFEC_PI*x[i+1]));
	s+=(x[m_numDim-1]-1)*(x[m_numDim-1]-1)*(1+sin(2*OFEC_PI*x[m_numDim-1])*sin(2*OFEC_PI*x[m_numDim-1]))+sin(3*OFEC_PI*x[0])*sin(3*OFEC_PI*x[0]);
	s=s*0.1;
	for(int i=0;i<m_numDim;i++)
		s+=u(x[i],5,100,4);

	obj[0]= s+m_bias;

}
double FPenalized_2::u(double x, double a, double k, double m)const{
	if(x>a) return k*pow(x-a,m);
	else if(x<-a) return k*pow(-x-a,m);
	else return 0;
}

