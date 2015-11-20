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
#include "FGriewank.h"


FGriewank::FGriewank(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
	if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank_noBounds_CEC05"))  {
          setSearchRange(0,600);
         for(int i=0;i<m_numDim; ++i)   m_searchRange[i].m_flag = false;
     }
	else if(IS_PROBLEM_NAME(m_id,"FUN_S_Griewank_Rosenbrock_F13_CEC05"))  {
          setSearchRange(-5, 5);
     }
     else setSearchRange(-600,600);

    initialize();
}
FGriewank::FGriewank(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
	if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank_noBounds_CEC05"))  {
          setSearchRange(0,600);
         for(int i=0;i<m_numDim; ++i)   m_searchRange[i].m_flag = false;
     }
	else if(IS_PROBLEM_NAME(m_id,"FUN_S_Griewank_Rosenbrock_F13_CEC05"))  {
          setSearchRange(-5, 5);
     }
     else setSearchRange(-600,600);

    initialize();
}

FGriewank::~FGriewank(){
    //dtor
}

void FGriewank::initialize(){
    setOriginalGlobalOpt();

	if(IS_PROBLEM_NAME(m_id,"FUN_Griewank")){
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Griewank_CEC08")){
		loadTranslation();
         setBias(-180);
	}else if(IS_PROBLEM_NAME(m_id,"FUN_R_Griewank")){
		setConditionNumber(2);
         loadRotation();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank")){
		setConditionNumber(3);
            loadTranslation();
            loadRotation();
		
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank_noBounds_CEC05")){
		    setConditionNumber(3);
            loadTranslation();
            loadRotation();
            setBias(-180);
	}else if(IS_PROBLEM_NAME(m_id,"FUN_S_Griewank_Rosenbrock_F13_CEC05")){
               loadTranslation();
               setBias(-130);	
	}else if(IS_PROBLEM_NAME(m_id,"FUN_Griewank_Rosenbrock_F13_CEC05")){
               setBias(0);	
	}else{
		 throw myException("Error: please check the problem ID@FGriewank::initialize");
	}

    setGlobalOpt();
    setAccuracy(1.0e-2);

}


void FGriewank::evaluate__(double const *x,vector<double>& obj){
     double result = 0;
	 if(IS_PROBLEM_NAME(m_id,"FUN_S_Griewank_Rosenbrock_F13_CEC05")  || IS_PROBLEM_NAME(m_id,"FUN_Griewank_Rosenbrock_F13_CEC05"))  {
          for(int i=0;i<m_numDim; ++i) {
               double result001 = 0;
               double result002 = 0;
               double x001 = x[i] + 1;
               double x002 = x[(i+1) % m_numDim] + 1;
               result001 += 100 * pow((x002 - x001 * x001), 2.0) + (x001 - 1) * (x001 - 1);
               result002 += result001 * result001 / 4000.0 - cos(result001 / sqrt((double)(i+1))) + 1;
               result += result002;
          }
          result += m_bias;
     }
     else {
          double s1=0,s2=1;
          for(int i=0;i<m_numDim;i++){
               s1+=x[i]*x[i]/4000.;
               s2*=cos(x[i]/sqrt((double) (i+1)));
          }
          result = s1-s2+1.+m_bias;
     }
     obj[0]= result;
}
