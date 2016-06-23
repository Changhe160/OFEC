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
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
/*******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
****************************************************************************************
*  LaTex:
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FModified_Rastrigin.h"


FModified_Rastrigin::FModified_Rastrigin(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1),m_k(2){
	v[param_numDim]=2;
	 setSearchRange(0.0,1.); // note
     initialize();
}
FModified_Rastrigin::FModified_Rastrigin(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1),m_k(2){
	 setSearchRange(0.0,1.0); // note
     initialize();
}

FModified_Rastrigin::~FModified_Rastrigin(){
    //dtor
}
void FModified_Rastrigin::initialize(){ // note
	m_k[0]=3;m_k[1]=4;
	setDisAccuracy(0.1);
	setAccuracy(1.e-5);
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(12); //production of m_ki
	ifstream in;
	stringstream ss;
	ss<<Global::g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/"<<m_name<<"_Opt_"<<m_numDim<<"D.txt";
	in.open(ss.str().c_str());
	if(in.fail()){
		throw myException("cannot open data file@FModified_Rastrigin::initialize()");
	}
	for(int i=0;i<12;++i){
		double x0,x1;
		in>>x0>>x1;
		m_globalOpt[i].data().m_x[0]=x0; m_globalOpt[i].data().m_x[1]=x1; m_globalOpt[i].data().m_obj[0]=2.0;
	}
	in.close();
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FModified_Rastrigin::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	for(int i=0;i<m_numDim;++i){
		s+=10+9*cos(2*OFEC_PI*m_k[i]*x[i]);
	}
	obj[0]= s+m_bias;  // note

}