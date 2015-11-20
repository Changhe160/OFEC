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
// Created: 11 May 2011
// Last modified:
#include "FSphere.h"

FSphere::FSphere(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
	setSearchRange(-100,100);
    initialize();

}
FSphere::FSphere(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-100,100);
    initialize();
}

FSphere::~FSphere(){
    //dtor
}
void FSphere::initialize(){
    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Sphere")||IS_PROBLEM_NAME(m_id,"FUN_Sphere_Noisy")||IS_PROBLEM_NAME(m_id,"FUN_Sphere_Noisy_CEC05")){}
	else if(IS_PROBLEM_NAME(m_id,"FUN_S_Sphere")){            loadTranslation();}
	else if(IS_PROBLEM_NAME(m_id,"FUN_R_Sphere")){ 
        setConditionNumber(2);
        loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Sphere")){
        setConditionNumber(2);
        loadTranslation();
        loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Sphere_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_S_Sphere_CEC08")){
        setBias(-450);
        loadTranslation();
	}else{
		throw myException("Error: please check the problem ID@FSphere::initialize");
		exit(0);
	}
    
    setGlobalOpt();

}


void FSphere::evaluate__(double const *x,vector<double>& obj){
    double fit=0;
    if(IS_PROBLEM_NAME(m_id,"FUN_Sphere_Noisy")){
        double noise;
        for(int i=0;i<m_numDim;i++){
			noise=0.01*Global::msp_global->mp_uniformPro->Next();
            fit+=(x[i]+noise)*(x[i]+noise);
        }
    }
    else if(IS_PROBLEM_NAME(m_id, "FUN_Sphere_Noisy_CEC05")){
        for(int i=0;i<m_numDim;i++)	fit+=x[i]*x[i];
		fit *= (1.0 + 0.1 * fabs(Global::msp_global->mp_normalPro->Next()));
    }
    else{
        for(int i=0;i<m_numDim;i++)	{ fit+=x[i]*x[i]; }
    }

    obj[0]= fit+m_bias;

}
