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
#include "FRosenbrock.h"

FRosenbrock::FRosenbrock(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-2.048,2.048);

     initialize();
}
FRosenbrock::FRosenbrock(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-2.048,2.048);

     initialize();
}

FRosenbrock::~FRosenbrock(){
    //dtor
}

void FRosenbrock::initialize(){
	vector<double> v(m_numDim,1);
    setOriginalGlobalOpt(0,&v);
	if(IS_PROBLEM_NAME(m_id,"FUN_Rosenbrock")){
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Rosenbrock")){
	loadTranslation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Rosenbrock_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_S_Rosenbrock_CEC08")){
		setBias(390);
            loadTranslation();
	}else{
		 throw myException("Error: please check the problem ID@FRosenbrock::initialize()");
	}
    setGlobalOpt(0,&v);
    setAccuracy(1.0e-2);
}

void FRosenbrock::evaluate__(double const *x,vector<double>& obj){

	double fitness=0;
	for(int i=0;i<m_numDim-1;i++){
		fitness+=100*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i])+(x[i]-1)*(x[i]-1);
	}

	obj[0]= fitness+m_bias;

}
