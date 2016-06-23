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
#include "FQuartic_Noisy.h"

FQuartic_Noisy::FQuartic_Noisy(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-1.28,1.28);

     initialize();
}
FQuartic_Noisy::FQuartic_Noisy(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-1.28,1.28);

     initialize();
}

FQuartic_Noisy::~FQuartic_Noisy(){
    //dtor
}
void FQuartic_Noisy::initialize(){
    setOriginalGlobalOpt();

    setGlobalOpt();
    setAccuracy(1.0e-2);

}

void FQuartic_Noisy::evaluate__(double const *x,vector<double>& obj){
   	double fitness=0;
	for(int i=0;i<m_numDim;i++){
		fitness+=(i+1)*pow(x[i],4);
	}
	fitness+=Global::msp_global->mp_uniformPro->Next();

    obj[0]= fitness+m_bias;

}
