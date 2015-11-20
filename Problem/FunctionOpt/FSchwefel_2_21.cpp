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
#include "FSchwefel_2_21.h"

FSchwefel_2_21::FSchwefel_2_21(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-100,100);

     initialize();
}
FSchwefel_2_21::FSchwefel_2_21(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-100,100);

     initialize();
}

FSchwefel_2_21::~FSchwefel_2_21(){
    //dtor
}
void FSchwefel_2_21::initialize(){
    setOriginalGlobalOpt();
    setGlobalOpt();

}

void FSchwefel_2_21::evaluate__(double const *x,vector<double>& obj){
   	double * y = new double[m_numDim];
	for(int i=0;i<m_numDim;i++) y[i]= fabs(x[i]);
	double max=*max_element(y,y+m_numDim);

	delete []y;
	y=0;
    obj[0]= max+m_bias;

}
