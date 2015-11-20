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
#include "FGear_Train.h"

FGear_Train::FGear_Train(ParamMap &v):Problem((v[param_proId]), 4,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 4,(v[param_proName]),1){
	 setSearchRange(12,60);
     initialize();

}

FGear_Train::FGear_Train(const int rId, const int rDim, string& rName):Problem(rId, 4, rName),\
	BenchmarkFunction(rId, 4, rName){

    setSearchRange(12,60);

     initialize();
}

FGear_Train::~FGear_Train(){
    //dtor
}


void FGear_Train::initialize(){
	vector<double> v(m_numDim,0);
	v[0]=15;v[1]=20;v[2]=57;v[3]=59;
	
    setOriginalGlobalOpt(0,&v);
    setGlobalOpt(0,&v);
}


void FGear_Train::evaluate__(double const *x,vector<double>& obj){

    int x1,x2,x3,x4;
	double s;
	x1=(int)x[0]; x2=(int)x[1]; x3=(int)x[2];x4=(int)x[3];
	s=1./6.931-x1*x2/(double)(x3*x4);

	obj[0]= s*s;

}
