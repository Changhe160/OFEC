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
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
/*******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
****************************************************************************************
*  LaTex:F_4(x)=\mathrm{e}^{-2\lg2*{{\frac{x-0.08}{0.854}}^2}}*{\sin{5\pi*{x^0.75-0.05}}}^6
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FBeasley.h"


FBeasley::FBeasley(ParamMap &v):Problem((v[param_proId]), 1,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 1,(v[param_proName]),1){
	
	v[param_numDim]=1;
	 setSearchRange(0,1.); // note
     initialize();
}
FBeasley::FBeasley(const int rId, const int rDimNumber, string& rName):Problem(rId, 1, rName,1),\
	BenchmarkFunction(rId, 1, rName,1){
	 setSearchRange(0,1.); // note
     initialize();
}

FBeasley::~FBeasley(){
    //dtor
}
void FBeasley::initialize(){ // note
	setOptType(MAX_OPT);  // note
	m_originalGlobalOpt.setNumOpts(5); //1 gopt+ 4 lopt
	setDisAccuracy(0.1);
	setAccuracy(1.e-6);
	m_originalGlobalOpt.setFlagLocTrue();

	m_originalGlobalOpt[0].data().m_x[0]=0.08;
	evaluate_(m_originalGlobalOpt[0].data(),false);
	m_originalGlobalOpt[1].data().m_x[0]=0.25;
	evaluate_(m_originalGlobalOpt[1].data(),false);
	m_originalGlobalOpt[2].data().m_x[0]=0.45;
	evaluate_(m_originalGlobalOpt[2].data(),false);
	m_originalGlobalOpt[3].data().m_x[0]=0.68;
	evaluate_(m_originalGlobalOpt[3].data(),false);
	m_originalGlobalOpt[4].data().m_x[0]=0.93;
	evaluate_(m_originalGlobalOpt[4].data(),false);

	m_globalOpt=m_originalGlobalOpt;
	addProTag(MMP);
}
void FBeasley::evaluate__(double const *x,vector<double>& obj){
	double s;
	s=exp(-2*log(2.)*((x[0]-0.08)/0.854)*((x[0]-0.08)/0.854))*pow(sin(5*OFEC_PI*(pow(x[0],0.75)-0.05)),6);
	obj[0]= s+m_bias;  // note

}