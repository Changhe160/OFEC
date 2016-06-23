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
// Created: 11 May 2011
// Last modified:
#include "FRastrigin.h"

FRastrigin::FRastrigin(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-5.12,5.12);

     initialize();
}
FRastrigin::FRastrigin(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-5.12,5.12);

     initialize();
}

FRastrigin::~FRastrigin(){
    //dtor
}
void FRastrigin::initialize(){
    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Rastrigin")||IS_PROBLEM_NAME(m_id,"FUN_Rastrigin_Noisy")){
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Rastrigin")){
		loadTranslation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Rastrigin")){
		setConditionNumber(2);
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Rastrigin")||IS_PROBLEM_NAME(m_id,"FUN_RS_Rastrigin_CEC05")){
		setConditionNumber(2);
            loadTranslation();
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Rastrigin_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_S_Rastrigin_CEC08")){
		 setBias(-330);
            loadTranslation();
	}else{
		 throw myException("Error: please check the problem ID@FRastrigin::initialize()");
	}

    setGlobalOpt();
    setAccuracy(1.0e-2);

}


void FRastrigin::evaluate__(double const *x,vector<double>& obj){

	double fit=0;
	if(IS_PROBLEM_NAME(m_id,"FUN_Rastrigin_Noisy")){
        double noise;
        for(int i=0;i<m_numDim;i++){
            noise=0.01*Global::msp_global->mp_uniformPro->Next();
            fit=fit+(x[i]+noise)*(x[i]+noise)-10.*cos(2*OFEC_PI*(x[i]+noise))+10.;
        }

	}else{
	    for(int i=0;i<m_numDim;i++)
		fit=fit+x[i]*x[i]-10.*cos(2*OFEC_PI*x[i])+10.;

	}

	obj[0]= fit+m_bias;


}
