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
#include "FParEst_FMSoundWaves.h"

FParEst_FMSoundWaves::FParEst_FMSoundWaves(ParamMap &v):Problem((v[param_proId]), 6,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 6,(v[param_proName]),1){
	 setSearchRange(-6.4,6.35);
     initialize();

}
FParEst_FMSoundWaves::FParEst_FMSoundWaves(const int rId,  const int rDim, string& rName):Problem(rId, 6, rName),\
	BenchmarkFunction(rId, 6, rName){

    setSearchRange(-6.4,6.35);

     initialize();
}

FParEst_FMSoundWaves::~FParEst_FMSoundWaves(){
    //dtor
}


void FParEst_FMSoundWaves::initialize(){
	vector<double> v(m_numDim);
	v[0]=1.0;v[1]=5.0;v[2]=1.5;v[3]=4.8; v[4]=2.0;v[5]=4.9;
	
    setOriginalGlobalOpt(0,&v);
    setGlobalOpt(0,&v);

}

void FParEst_FMSoundWaves::evaluate__(double const *x,vector<double>& obj){

	double theta=2*OFEC_PI/100.;
	double s=0,t;
	for(int i=0;i<100;i++){
		t=x[0]*sin(x[1]*i*theta+x[2]*sin(x[3]*i*theta+x[4]*sin(x[5]*i*theta)))-sin(5.*i*theta+1.5*sin(4.8*i*theta+2.0*sin(4.9*i*theta)));
		s+=t*t;
	}
	obj[0]= s;

}
