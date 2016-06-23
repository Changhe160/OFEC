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
#include "FWaves.h"

FWaves::FWaves(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1){
	vector<double> lower, upper;
	 v[param_numDim]=2;
     upper.push_back (1.2);
	 upper.push_back (1.2);
	 lower.push_back (-0.9);
	 lower.push_back (-1.2);
	 setSearchRange(lower,upper);
     initialize();
}
FWaves::FWaves(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1){
	vector<double> lower, upper;
	 upper.push_back (1.2);
	 upper.push_back (1.2);
	 lower.push_back (-0.9);
	 lower.push_back (-1.2);
	 setSearchRange(lower,upper);
     initialize();
}

FWaves::~FWaves(){
    //dtor
}
void FWaves::initialize(){
	setOptType(MAX_OPT); 
	setDisAccuracy(0.15);
	setAccuracy(1.e-3);
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(10);	//one global optimum+9 local optimum

	ifstream in;
	stringstream ss;
	ss<<Global::g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/"<<m_name<<"_Opt_"<<m_numDim<<"Dim.txt";
	in.open(ss.str().c_str());
	if(!in)		throw myException("cannot open data file@FShubert::initialize()");

	for(int i=0;i<10;++i){
		double x0,x1,obj;
		in>>x0>>x1>>obj;
		m_globalOpt[i].data().m_x[0]=x0; m_globalOpt[i].data().m_x[1]=x1; m_globalOpt[i].data().m_obj[0]=obj;
	}
	in.close();
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FWaves::evaluate__(double const *x,vector<double>& obj){
	double s,t;
	t=x[1]*x[1];
	s=pow(0.3*x[0],3)+3.5*x[0]*pow(x[1],3)-4.7*cos(3*x[0]-(2+x[0])*t)*sin(2.5*x[0]*OFEC_PI);
	obj[0]= s+m_bias;
}