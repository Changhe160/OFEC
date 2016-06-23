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
#include "FMAX_global5.h"


FMAX_global5::FMAX_global5(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	 v[param_numDim]=2;
    setSearchRange(-6,6);

     initialize();
}
FMAX_global5::FMAX_global5(const int rId,  const int rDim, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){

    setSearchRange(-6,6);

     initialize();
}

FMAX_global5::~FMAX_global5(){
    //dtor
}


void FMAX_global5::initialize(){
	m_OptMode[0]=MAX_OPT;
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(5); //4 gopt + 1 lopt
	setDisAccuracy(0.1);
	setAccuracy(1.e-5);

	CodeVReal x(2,1);
	x.m_x[0]=3.58443; x.m_x[1]=-1.84813; x.m_obj[0]=200.;
	m_globalOpt[0].data()=x;
	x.m_x[0]=3.; x.m_x[1]=2.; x.m_obj[0]=200.;
	m_globalOpt[1].data()=x;
	x.m_x[0]=-2.80512; x.m_x[1]=3.13131; x.m_obj[0]=200.;
	m_globalOpt[2].data()=x;
	x.m_x[0]=-3.77931; x.m_x[1]=-3.28319; x.m_obj[0]=200.;
	m_globalOpt[3].data()=x;
	x.m_x[0]=3.49157; x.m_x[1]=6.; x.m_obj[0]=-907.413;
	m_globalOpt[4].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}

void FMAX_global5::evaluate__(double const *x,vector<double>& obj){

   obj[0]= 200-pow(x[0]*x[0]+x[1]-11,2.)-pow(x[0]+x[1]*x[1]-7,2.)+m_bias;

}

