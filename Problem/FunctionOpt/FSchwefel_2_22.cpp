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
#include "FSchwefel_2_22.h"

FSchwefel_2_22::FSchwefel_2_22(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-10,10);

     initialize();
}
FSchwefel_2_22::FSchwefel_2_22(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-10,10);

     initialize();
}

FSchwefel_2_22::~FSchwefel_2_22(){
    //dtor
}

void FSchwefel_2_22::initialize(){

    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_2_22")||IS_PROBLEM_NAME(m_id,"FUN_Schwefel_2_22_Noisy")){
	
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_2_22")){
		 loadTranslation();
	
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Schwefel_2_22")){
		 setConditionNumber(2);
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Schwefel_2_22")){
		setConditionNumber(2);
            loadTranslation();
            loadRotation();
	}else {
		 throw myException("Error: please check the problem ID@FSchwefel_2_22::initialize");
	}
    setGlobalOpt();

}

void FSchwefel_2_22::evaluate__(double const *x,vector<double>& obj){

    double s1=0,s2=1.;
    if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_2_22_Noisy")){
        double noise;
         for(int i=0;i<m_numDim;i++){
            noise=0.01*Global::msp_global->mp_uniformPro->Next();
            s1+=fabs(x[i]+noise);
            s2*=fabs(x[i]+noise);
        }

    }else{
        for(int i=0;i<m_numDim;i++){
            s1+=fabs(x[i]);
            s2*=fabs(x[i]);
        }
    }

	obj[0]= s1+s2+m_bias;


}
