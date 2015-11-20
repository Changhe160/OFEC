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
#include "Algorithm.h"
#include "../Global/global.h"
#include "../Problem/problem.h"
#ifdef OFEC_DEMON
extern bool g_algTermination;
#endif
Algorithm::~Algorithm()
{
    //dtor
}
Algorithm::Algorithm(const int rID, string rName):m_algID(rID),m_name(rName),\
	m_numDim((Global::msp_global->mp_problem.get()!=nullptr)?GET_NUM_DIM:0),m_isTerminated(false){

}
Algorithm & Algorithm::operator=(const Algorithm& rhs){
	if(m_algID!=rhs.m_algID) {
		throw myException("Algorithm assignment@Algorithm::operator=(const Algorithm& rhs)");
		return *this;
	}
	m_name=rhs.m_name;
	m_numDim=rhs.m_numDim;
	m_isTerminated=rhs.m_isTerminated;
	m_algPar.str(rhs.m_algPar.str());
	return *this;
}

bool Algorithm::ifTerminating(){
	#ifdef OFEC_CONSOLE
	if(Global::g_arg.find(param_maxEvals)!= Global::g_arg.end())
	if(Global::msp_global->mp_problem->getEvaluations()>=(int)Global::g_arg[param_maxEvals]) return true;	
	#endif

	#ifdef OFEC_DEMON
	if(g_algTermination) return true;
	#endif

	return false;
}
bool Algorithm::ifTerminated(){
	return m_isTerminated;
}
ReturnFlag Algorithm::run(){
	ReturnFlag rf=run_();
	m_isTerminated=true;
	return rf;
};

ReturnFlag Algorithm::run_(){
	return Return_Terminate;
}