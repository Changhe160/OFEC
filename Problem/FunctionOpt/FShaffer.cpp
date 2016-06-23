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
#include "FShaffer.h"

FShaffer::FShaffer(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]),2,(v[param_proName]),1){
	v[param_numDim]=2;
	 setSearchRange(-15,15);
     initialize();
}
FShaffer::FShaffer(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	 setSearchRange(-15,15);
     initialize();
}

FShaffer::~FShaffer(){
    //dtor
}
void FShaffer::initialize(){
	setAccuracy(1.e-6);
	setDisAccuracy(0.1);
	setOptType(MAX_OPT);
	 m_globalOpt.setFlagLocTrue();
	 m_originalGlobalOpt.setFlagLocTrue();
	 m_globalOpt.setNumOpts(1);
	 m_originalGlobalOpt.setNumOpts(1);

	 CodeVReal x(m_numDim,1);
	 x.m_x[0]=0.0;x.m_x[1]=0.0; x.m_obj[0]=0.9999;
	 m_globalOpt[0].data()=m_originalGlobalOpt[0].data()=x;
	 setObjSet();
}
void FShaffer::evaluate__(double const *x,vector<double>& obj){
	double s,t=x[0]*x[0]+x[1]*x[1];
	s=0.5+(0.5-pow(sin(sqrt(0.0001+t)),2))/pow(1+0.001*t*t,2);
	obj[0]= s+m_bias;

}