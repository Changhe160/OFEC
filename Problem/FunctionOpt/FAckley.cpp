
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
#include "FAckley.h"


FAckley::FAckley(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-32.768,32.768);
     initialize();
}
FAckley::FAckley(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
    setSearchRange(-32.768,32.768);
     initialize();
}

FAckley::~FAckley(){
    //dtor
}

void FAckley::initialize(){

    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Ackley")||IS_PROBLEM_NAME(m_id,"FUN_Ackley_Noisy")){	}
	else if(IS_PROBLEM_NAME(m_id,"FUN_S_Ackley")){
		loadTranslation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Ackley")){
		 setConditionNumber(2);
         loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Ackley")){
		setConditionNumber(100);
        loadTranslation();
        loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Ackley_Bound_CEC05")){
		setConditionNumber(100);
        loadTranslation();
        loadRotation();
        setBias(-140);

	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Ackley_CEC08")){
		setBias(-140);
            loadTranslation();
	}else{
		throw myException("Error: please check the problem ID@FAckley::initialize()");
	}

    setGlobalOpt();
}
 
void FAckley::evaluate__(double const *x,vector<double>& obj){

	double s1=0,s2=0;

	if(IS_PROBLEM_NAME(m_id,"FUN_Ackley_Noisy")){

        double noise;
        for(int i=0;i<m_numDim;i++){
            noise=0.01*Global::msp_global->mp_uniformPro->Next();
            s1+=(x[i]+noise)*(x[i]+noise);
            s2+=cos(2*OFEC_PI*(x[i]+noise));
        }
	}else{
        for(int i=0;i<m_numDim;i++){
		s1+=x[i]*x[i];
		s2+=cos(2*OFEC_PI*x[i]);
        }
	}

	obj[0]= -20*exp(-0.2*sqrt(s1/m_numDim))-exp(s2/m_numDim)+20+OFEC_E+m_bias;

}
