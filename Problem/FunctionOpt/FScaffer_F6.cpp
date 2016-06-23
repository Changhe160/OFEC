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
#include "FScaffer_F6.h"

FScaffer_F6::FScaffer_F6(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    setSearchRange(-100,100);
     initialize();
}
FScaffer_F6::FScaffer_F6(const int rId, const int rDim, string& rName):Problem(rId, rDim, rName,1),\
	BenchmarkFunction(rId, rDim, rName,1){

    setSearchRange(-100,100);
     initialize();
}
FScaffer_F6::~FScaffer_F6(){
    //dtor
}

void FScaffer_F6::initialize(){
     setOriginalGlobalOpt();
	 if(IS_PROBLEM_NAME(m_id,"FUN_RS_Expanded_Scaffer_F6_CEC05")){
		setConditionNumber(3);
            loadTranslation();
            loadRotation();
            setBias(-300);
	 }else{
		setBias(0);
	 }
    
    setGlobalOpt();
    setAccuracy(1.0e-2);
}

void FScaffer_F6::evaluate__(double const *x,vector<double>& obj){
     double fitness = 0;
     if(  IS_PROBLEM_NAME(m_id,"FUN_RS_Expanded_Scaffer_F6_CEC05") || IS_PROBLEM_NAME(m_id,"FUN_Expanded_Scaffer_F6_CEC05")) {
          for(int i=0;i<m_numDim; ++i) {
               double result001 = 0;
               double x001 = x[i];
               double x002 = x[(i+1) % m_numDim];
               result001 = 0.5 + (pow(sin(sqrt(x001*x001+x002*x002)), 2.0)-0.5)/ pow((1+0.001*(x001*x001+x002*x002)), 2.0);
               fitness += result001;
          }
     }
     else if(IS_PROBLEM_NAME(m_id,"FUN_Noncont_Expanded_Scaffer_F6_CEC05") ) {
          double *y=new double[m_numDim];
          for(int i=0;i<m_numDim;i++){
               if(fabs(x[i])<0.5) y[i]=x[i];
               else {
                    double a,b=2*x[i];
                    if(b <= 0 && b-(int)b<0.5) a=(int)b;
                    else if(b-(int)b<0.5) a=(int)b;
                    else a = (int)b + 1;
                    y[i]=a/2;
               }
               double result001 = 0;
               double x001 = y[i];
               double x002 = y[(i+1) % m_numDim];
               result001 = 0.5 + (pow(sin(sqrt(x001*x001+x002*x002)), 2.0)-0.5)/ pow((1+0.001*(x001*x001+x002*x002)), 2.0);
               fitness += result001;
          }
          delete[] y;
     }
	else {
          fitness=0.5+(sin(sqrt(x[0]*x[0]+x[1]*x[1]))*sin(sqrt(x[0]*x[0]+x[1]*x[1]))-0.5)/((1+0.001*(x[0]*x[0]+x[1]*x[1]))*(1+0.001*(x[0]*x[0]+x[1]*x[1])));
	}

	 obj[0]= fitness+m_bias;

}
