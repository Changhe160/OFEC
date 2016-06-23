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
******************************************************************************************
*  Paper; Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FKeane_Bump.h"


FKeane_Bump::FKeane_Bump(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
     
	 setSearchRange(0,10);
     initialize();
}
FKeane_Bump::FKeane_Bump(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
	 
	 setSearchRange(0,10);
     initialize();
}

FKeane_Bump::~FKeane_Bump(){
    //dtor
}
void FKeane_Bump::initialize(){
	setOptType(MAX_OPT);
	setDisAccuracy(0.5);
	setAccuracy(1.e-4);
	m_globalOpt.flagGloObj()=false;
	m_globalOpt.flagLoc()=false;
	m_globalOpt.setNumOpts(1);
	m_originalGlobalOpt=m_globalOpt;
	setObjSet();
}
bool FKeane_Bump::isValid(const VirtualEncoding &x_){
	const CodeVReal&x = dynamic_cast<const CodeVReal&>(x_);
	 double m=1,n=0;
	 int i;
	 for(i=0;i<m_numDim;i++){
	    m*=x[i];
		n+=x[i];
	 }
	 if(m>0.75&&n<15*m_numDim/2) return true;
	 else return false;
};
void FKeane_Bump::evaluate__(double const *x,vector<double>& obj){
	double s,a=0,b=1,c=0;
	int i;
	for(i=0;i<m_numDim;i++){
		a+=pow(cos(x[i]),4);
		b*=cos(x[i])*cos(x[i]);
		c+=i*x[i]*x[i];
	}
	s=abs(a-2*b)/sqrt(c);	
	obj[0]= s;

}