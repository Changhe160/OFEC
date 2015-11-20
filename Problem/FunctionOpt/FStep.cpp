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
#include "FStep.h"

FStep::FStep(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-100,100);

     initialize();
}
FStep::FStep(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-100,100);

     initialize();
}

FStep::~FStep(){
    //dtor
}
void FStep::initialize(){
    setOriginalGlobalOpt();

    setGlobalOpt();

}

void FStep::evaluate__(double const *x,vector<double>& obj){
    double fitness=0;
	for(int i=0;i<m_numDim;i++){
		fitness+=fabs((double)int(x[i]+0.5)*int(x[i]+0.5));
	}

    obj[0]= fitness+m_bias;

}

