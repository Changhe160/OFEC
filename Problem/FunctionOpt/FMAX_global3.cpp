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
#include "FMAX_global3.h"


FMAX_global3::FMAX_global3(ParamMap &v):Problem((v[param_proId]), 1,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 1,(v[param_proName]),1){
	v[param_numDim]=1;
	setSearchRange(0,1);
	initialize();
}
FMAX_global3::FMAX_global3(const int rId,  const int rDim,  string& rName):Problem(rId, 1, rName,1),\
	BenchmarkFunction(rId, 1, rName,1){
	setSearchRange(0,1);
	initialize();
}

FMAX_global3::~FMAX_global3(){
    //dtor
}

void FMAX_global3::initialize(){
    m_OptMode[0]=MAX_OPT;
    m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(5); //3 gopt + 2 lopt
	setDisAccuracy(0.1);
	setAccuracy(1.e-5);

	CodeVReal x(1,1);
	x.m_x[0]=0.0; x.m_obj[0]=0.125;
	m_globalOpt[0].data()=x;
	x.m_x[0]=0.933895; x.m_obj[0]=1.;
	m_globalOpt[1].data()=x;
	x.m_x[0]=0.246655; x.m_obj[0]=1.;
	m_globalOpt[2].data()=x;
	x.m_x[0]=0.694574; x.m_obj[0]=0.91642;
	m_globalOpt[3].data()=x;
	x.m_x[0]=0.450627; x.m_obj[0]=1.;
	m_globalOpt[4].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
}


void FMAX_global3::evaluate__(double const *x,vector<double>& obj){

	obj[0]= pow(sin(5*OFEC_PI*(pow(x[0],3./4.)-0.05)),6.)+m_bias;

}

