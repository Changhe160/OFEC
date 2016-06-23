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
#include "FSchwefel_1_2.h"

FSchwefel_1_2::FSchwefel_1_2(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-100,100);

     initialize();
}
FSchwefel_1_2::FSchwefel_1_2(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){

    setSearchRange(-100,100);

     initialize();
}

FSchwefel_1_2::~FSchwefel_1_2(){
    //dtor
}

void FSchwefel_1_2::initialize(){

    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_1_2")||IS_PROBLEM_NAME(m_id,"FUN_Schwefel_1_2_Noisy")){
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2")||IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2_Noisy")){
		loadTranslation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Schwefel_1_2")){
		setConditionNumber(2);
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Schwefel_1_2")){
		 setConditionNumber(2);
            loadTranslation();
            loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2_Noisy_CEC05")){
		 setBias(-450);
            loadTranslation();
	}else{
		throw myException("Error: please check the problem ID@FSchwefel_1_2::initialize");
	}

    setGlobalOpt();
    setAccuracy(1.0e-2);
}


void FSchwefel_1_2::evaluate__(double const *x,vector<double>& obj){

    double s1=0,s2=0;

    if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_1_2_Noisy") ||IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2_Noisy")){
        double noise;
        for(int i=0;i<m_numDim;i++){
            for(int j=0;j<=i;j++){
                noise=0.01*Global::msp_global->mp_uniformPro->Next();
                s1+=(x[j]+noise);
            }
            s2+=s1*s1;
            s1=0;
        }
    }
    else if(IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel_1_2_Noisy_CEC05") ) {
          for(int i=0;i<m_numDim;i++){
               for(int j=0;j<=i;j++)
                    s1+=x[j];
               s2+=s1*s1;
               s1=0;
          }
		  s2 *= (1.0 + 0.4 * fabs(Global::msp_global->mp_normalPro->Next()));
    }
    else{
        for(int i=0;i<m_numDim;i++){
		for(int j=0;j<=i;j++)
			s1+=x[j];
		s2+=s1*s1;
		s1=0;
        }
    }

 obj[0]= s2+m_bias;

}
