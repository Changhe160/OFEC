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
#include "FSchwefel.h"


FSchwefel::FSchwefel(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1)
{

    setSearchRange(-500,500);

    initialize();
}
FSchwefel::FSchwefel(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-500,500);

    initialize();
}

FSchwefel::~FSchwefel()
{
    //dtor
}


void FSchwefel::initialize()
{
	vector<double> v(m_numDim,420.9687);
    setOriginalGlobalOpt(0,&v);

	if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel")||IS_PROBLEM_NAME(m_id,"FUN_Schwefel_Noisy")){
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel")){
		loadTranslation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Schwefel")){
		setConditionNumber(2);
        loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Schwefel")){
		 setConditionNumber(2);
        loadTranslation();
        loadRotation();
	}else{
		throw myException("Error: please check the problem ID@FSchwefel::initialize");
	}
	
    setGlobalOpt(0,&v);
    setAccuracy(1.0e-2);

}

void FSchwefel::evaluate__(double const *x,vector<double>& obj)
{

    double fitness=0;
    for(int i=0;i<m_numDim; i++)
    {
        fitness+=-x[i]*sin(sqrt(fabs(x[i])));
    }

    obj[0]= fitness+m_bias;

}
