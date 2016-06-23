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
#include "FNoncont_Rastrigin.h"

FNoncont_Rastrigin::FNoncont_Rastrigin(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
    setSearchRange(-5.12,5.12);
    initialize();
}
FNoncont_Rastrigin::FNoncont_Rastrigin(const int rId, const int rDim, string& rName):Problem(rId, rDim, rName,1),\
	BenchmarkFunction(rId, rDim, rName,1){
    setSearchRange(-5.12,5.12);
    initialize();
}

FNoncont_Rastrigin::~FNoncont_Rastrigin(){
    //dtor
}

void FNoncont_Rastrigin::initialize(){
    setOriginalGlobalOpt();
    setGlobalOpt();
}

void FNoncont_Rastrigin::evaluate__(double const *x,vector<double>& obj){

    double fit=0;
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
	}
	for(int i=0;i<m_numDim;i++){
		fit=fit+y[i]*y[i]-10.*cos(2*OFEC_PI*y[i])+10.;
	}
	delete [] y;
	y=0;

    obj[0]= fit+m_bias;

}
