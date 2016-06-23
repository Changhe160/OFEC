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
#include "FMAX_global1.h"


FMAX_global1::FMAX_global1(ParamMap &v):Problem((v[param_proId]), 1,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 1,(v[param_proName]),1){
	v[param_numDim]=1;
    setSearchRange(0,1);
    initialize();
}
FMAX_global1::FMAX_global1(const int rId,  const int rDim, string& rName):Problem(rId, 1, rName,1),\
	BenchmarkFunction(rId, 1, rName,1){
    setSearchRange(0,1);
    initialize();
}

FMAX_global1::~FMAX_global1(){
    //dtor
}

void FMAX_global1::initialize(){
    m_OptMode[0]=MAX_OPT;
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(5); //5 gopt
	setDisAccuracy(0.1);
	setAccuracy(1.e-5);
	
	CodeVReal x(1,1);
	x.m_x[0]=0.5; x.m_obj[0]=1.;
	m_globalOpt[0].data()=x;
	x.m_x[0]=0.1;
	m_globalOpt[1].data()=x;
	x.m_x[0]=0.3;
	m_globalOpt[2].data()=x;
	x.m_x[0]=0.7;
	m_globalOpt[3].data()=x;
	x.m_x[0]=0.9;
	m_globalOpt[4].data()=x;
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}

void FMAX_global1::evaluate__(double const *x,vector<double>& obj){

	obj[0]=  pow(sin(5*OFEC_PI*x[0]),6.)+m_bias;

}
