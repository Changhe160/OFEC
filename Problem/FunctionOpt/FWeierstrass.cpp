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
#include "FWeierstrass.h"

FWeierstrass::FWeierstrass(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-0.5,0.5);

     initialize();
}
FWeierstrass::FWeierstrass(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-0.5,0.5);

     initialize();
}

FWeierstrass::~FWeierstrass(){
    //dtor
}


void FWeierstrass::initialize(){

    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Weierstrass")){
	
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Weierstrass")){
		setConditionNumber(2);
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Weierstrass")){
		setConditionNumber(5);
            loadTranslation();
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Weierstrass_CEC05")){
		setConditionNumber(5);
            loadTranslation();
            loadRotation();
             setBias(90);
	}else {
		throw myException("Error: please check the problem ID@FWeierstrass::initialize");
	}

    setGlobalOpt();
    setAccuracy(1.0e-2);

}


void FWeierstrass::evaluate__(double const *x,vector<double>& obj){

	double a=0.5,b=3;
	int kmax=20;
	double fit=0,s=0;
	for(int i=0;i<m_numDim;i++)
		for(int k=0;k<=kmax;k++)
			fit+=pow(a,k)*cos(2*OFEC_PI*pow(b,k)*(x[i]+0.5));

	for(int k=0;k<=kmax;k++)
			s+=pow(a,k)*cos(2*OFEC_PI*pow(b,k)*0.5);
	s=s*m_numDim;
	obj[0]= fit-s+m_bias;

}
