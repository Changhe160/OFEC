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
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FSzu.h"
//* 
//* F(vec3{X})=\sum_{i=1}^{D}{-x_i^2}
//*

FSzu::FSzu(ParamMap &v):Problem((v[param_proId]),(v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]),(v[param_numDim]),(v[param_proName]),1){
	 setSearchRange(-5.0,5.0);
     initialize();
}
FSzu::FSzu(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
	 setSearchRange(-5.,5.);
     initialize();
}

FSzu::~FSzu(){
    //dtor
}
void FSzu::initialize(){
	setOptType(MIN_OPT);

	 m_originalGlobalOpt.setNumOpts(1);
	 vector<vector<double>> gobj;

	 if(m_numDim==2){
		 gobj.push_back(vector<double>(1,-156.66));	 
	 }else if(m_numDim==3) {
		 gobj.push_back(vector<double>(1,-235.0));
	 }else if(m_numDim==4) {
		 gobj.push_back(vector<double>(1,-313.33));
	 }else if(m_numDim==5) {
		 gobj.push_back(vector<double>(1,-391.66));
	 }else if(m_numDim==6) {
		 gobj.push_back(vector<double>(1,-469.99));
	 }else if(m_numDim==7) {
		 gobj.push_back(vector<double>(1,-548.33));	
	 }else if(m_numDim==8) {
		 gobj.push_back(vector<double>(1,-626.66));	
	 }else if(m_numDim==9) {
		 gobj.push_back(vector<double>(1,-704.99));		
	 }

	 if(m_numDim>=2&&m_numDim<=9){
		 m_originalGlobalOpt.setGloObj(gobj);
		 m_originalGlobalOpt.flagGloObj()=true;
		 m_originalGlobalOpt.flagLoc()=false;
	 }
	 else{
		m_originalGlobalOpt.flagGloObj()=false;
		m_originalGlobalOpt.flagLoc()=false;
	 } 
	 m_globalOpt=m_originalGlobalOpt;
	 setDisAccuracy(0.1);
	 setAccuracy(1.e-2);
}
void FSzu::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	int i;
	for(i=0;i<m_numDim;++i){
		s+=pow(x[i],4)-16*x[i]*x[i]+5*x[i];
	}
	obj[0]= s;

}