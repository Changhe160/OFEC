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
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FHimmenblau.h"


FHimmenblau::FHimmenblau(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	 v[param_numDim]=2;
	 setSearchRange(-6,6);
     initialize();
}
FHimmenblau::FHimmenblau(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	 setSearchRange(-6,6);
     initialize();
}

FHimmenblau::~FHimmenblau(){
    //dtor
}
void FHimmenblau::initialize(){
	 setDisAccuracy(0.5);
	 setAccuracy(1.e-4);
	 m_globalOpt.setFlagLocTrue();
	 m_globalOpt.setNumOpts(4); // 1 gopt+3 lopt

	 CodeVReal x(m_numDim,1);
	 x.m_x[0]=3.0;x.m_x[1]=2.0;
	 evaluate_(x,false);
	 m_globalOpt[0].data()=x;
	 x.m_x[0]=3.58149;x.m_x[1]=-1.8208 ; 
	 evaluate_(x,false);
	 m_globalOpt[1].data()=x;
	 x.m_x[0]=-2.78706 ;x.m_x[1]=3.1282 ; 
	 evaluate_(x,false);
	 m_globalOpt[2].data()=x;
	 x.m_x[0]=-3.76343;x.m_x[1]=-3.26605; 
	 evaluate_(x,false);
	 m_globalOpt[3].data()=x;
	 
	 m_originalGlobalOpt=m_globalOpt;
	 addProTag(MMP);

	 setObjSet();
}
void FHimmenblau::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	double t0=(x[0]*x[0]+x[1]-11),t1=(x[1]*x[1]+x[0]-7);
	s=t0*t0+t1*t1+0.1*(pow(x[0]-3,2)+pow(x[1]-2,2));
	obj[0]= s+m_bias;
}