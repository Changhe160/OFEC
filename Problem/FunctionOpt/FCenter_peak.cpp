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
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#include "FCenter_peak.h"

FCenter_peak::FCenter_peak(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	 v[param_numDim]=2;
	 setSearchRange(-2,2);
     initialize();
}
FCenter_peak::FCenter_peak(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	
	 setSearchRange(-2,2);
     initialize();
}

FCenter_peak::~FCenter_peak(){
    //dtor
}
void FCenter_peak::initialize(){
	 setOptType(MAX_OPT); 
	 m_globalOpt.setFlagLocTrue();
	 m_globalOpt.setNumOpts(5); //1 gopt + 4 lopt
	 setDisAccuracy(0.2);
	 setAccuracy(1.e-5);

	ifstream in;
	stringstream ss;
	ss<<Global::g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/"<<m_name<<"_Opt_"<<m_numDim<<"Dim.txt";
	in.open(ss.str().c_str());
	if(!in)		throw myException("cannot open data file@FShubert::initialize()");

	for(int i=0;i<5;++i){
		double x0,x1,obj;
		in>>x0>>x1>>obj;
		m_globalOpt[i].data().m_x[0]=x0; m_globalOpt[i].data().m_x[1]=x1; m_globalOpt[i].data().m_obj[0]=obj;
	}
	in.close();
	
	 m_originalGlobalOpt=m_globalOpt;
	 addProTag(MMP);
	 setObjSet();
}
void FCenter_peak::evaluate__(double const *x,vector<double>& obj){
	double s;
	s=3*sin(0.5*OFEC_PI*x[0]+0.5*OFEC_PI)*(2-sqrt(x[0]*x[0]+x[1]*x[1]))/4;
	obj[0]= s+m_bias;

}