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
#include "FElliptic.h"

FElliptic::FElliptic(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-100.,100.);

     initialize();
}
FElliptic::FElliptic(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-100.,100.);

     initialize();
}

FElliptic::~FElliptic(){
    //dtor
}

void FElliptic::initialize(){

    setOriginalGlobalOpt();
	if(IS_PROBLEM_NAME(m_id,"FUN_Elliptic")){

	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Elliptic_CEC05")){
            loadTranslation();
            loadRotation();
            setBias(-450);	
	}else{
		throw myException("Error: please check the problem ID@FElliptic::initialize");
	}

    setGlobalOpt();

}



void FElliptic::evaluate__(double const *x,vector<double>& obj){

	double s=0;
	for(int i=0;i< m_numDim;i++){
		s+=pow(1e6,i/(m_numDim-1.))*x[i]*x[i];
	}
    obj[0]= s+m_bias;
}
