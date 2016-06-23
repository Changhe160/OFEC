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
#include "FMAX_global4.h"


FMAX_global4::FMAX_global4(ParamMap &v):Problem((v[param_proId]), 1,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 1,(v[param_proName]),1){
	v[param_numDim]=1;
	setSearchRange(0,1);
	initialize();
}
FMAX_global4::FMAX_global4(const int rId,  const int rDim, string& rName):Problem(rId, 1, rName,1),\
	BenchmarkFunction(rId, 1, rName,1){
	setSearchRange(0,1);
	initialize();
}

FMAX_global4::~FMAX_global4(){
    //dtor
}


void FMAX_global4::initialize(){
    m_OptMode[0]=MAX_OPT;
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(5); //1 gopt + 4 lopt
	setDisAccuracy(0.1);
	setAccuracy(1.e-5);

	CodeVReal x(1,1);
	x.m_x[0]=0.0796998; x.m_obj[0]=1.;
	m_globalOpt[0].data()=x;
	x.m_x[0]=0.930153; x.m_obj[0]=0.25161;
	m_globalOpt[1].data()=x;
	x.m_x[0]=0.246279; x.m_obj[0]=0.94868;
	m_globalOpt[2].data()=x;
	x.m_x[0]=0.679166; x.m_obj[0]=0.50411;
	m_globalOpt[3].data()=x;
	x.m_x[0]=0.449496; x.m_obj[0]=0.77081;
	m_globalOpt[4].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}

void FMAX_global4::evaluate__(double const *x,vector<double>& obj){

    double s;
	s=pow((x[0]-0.08)/.854,2.);
	s=-s*2.*log(2.);

	obj[0]= exp(s)*pow(sin(5*OFEC_PI*(pow(x[0],3./4.)-0.05)),6.)+m_bias;

}
