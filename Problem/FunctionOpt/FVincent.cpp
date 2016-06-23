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
/*******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
****************************************************************************************
*  LaTex:F(x)=\frac{1}{D}\sum^D_{i=1}{sin(10\log(x_i))}
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FVincent.h"

FVincent::FVincent(ParamMap &v):Problem((v[param_proId]), v[param_numDim],(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), v[param_numDim],(v[param_proName]),1){
	
	 setSearchRange(0.25,10.); // note
     initialize();
}
FVincent::FVincent(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
	 setSearchRange(0.25,10.); // note
     initialize();
}

FVincent::~FVincent(){
    //dtor
}
void FVincent::initialize(){ // note
	setOptType(MAX_OPT);  // note
	setDisAccuracy(0.2);
	setAccuracy(1.e-4);

	m_globalOpt.flagLoc()=(false);
	m_globalOpt.flagGloObj()=(true);
	m_globalOpt.setNumOpts(pow(6,m_numDim));
	m_globalOpt.setGloObj(vector<vector<double>>(pow(6,m_numDim),vector<double>(m_numObj,1)));
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FVincent::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	for(int i=0;i<m_numDim;++i){
		s+=sin(10*log(x[i]));
	}
	s/=m_numDim;
	obj[0]= s+m_bias;  // note

}