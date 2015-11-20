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
// Created: 11 May 2011
// Last modified:

#include "problem.h"
#include "../Global/global.h"
Problem::Problem(const int rId, const int rDimNumber, string rName, const int numObj):m_id(rId),\
m_numDim(rDimNumber), m_OptMode(numObj,MIN_OPT), m_evals(0), m_numObj(numObj), m_name(rName), \
	m_popInitialMode(POP_INIT_UNIFORM),m_validationMode(VALIDATION_REMAP),m_accuracy(0),m_sameType(true),\
	m_tag(1,SOP),m_tevals(0),m_cevals(0){
    
	int mode=Global::g_arg[param_solutionValidationMode];	
	m_validationMode=static_cast<SolutionValidation>(mode);
	
	mode=Global::g_arg[param_populationInitialMethod];
	m_popInitialMode=static_cast<PopInitMethod>(mode);

    m_proPar<<"Dimension:"<<m_numDim<<"; Number of Objectives: "<<m_numObj<<"; ";

}

Problem& Problem::operator=(const Problem & rhs){
    if(this== &rhs) return *this;

	m_OptMode=rhs.m_OptMode;
    m_evals=rhs.m_evals;
	m_numObj=rhs.m_numObj;
	m_name=rhs.m_name;   
	m_validationMode=rhs.m_validationMode;
	m_popInitialMode=rhs.m_popInitialMode;
	m_accuracy=rhs.m_accuracy;	
	m_sameType=rhs.m_sameType;

	m_tag=rhs.m_tag;
	m_tevals=rhs.m_tevals;
	m_cevals=rhs.m_cevals;

	m_numDim = rhs.m_numDim;

    return *this;
}
void Problem::parameterSetting(Problem * rhs){
	m_accuracy=rhs->m_accuracy;
	m_OptMode=rhs->m_OptMode;
    m_evals=rhs->m_evals;
	m_sameType=rhs->m_sameType;
	m_validationMode=rhs->m_validationMode;
	m_popInitialMode=rhs->m_popInitialMode;
	m_tag=rhs->m_tag;
	m_tevals=rhs->m_tevals;
	m_cevals=rhs->m_cevals;
}

void Problem::setAccuracy(double rAcc){
    m_accuracy=rAcc;
}
void Problem::setObjNumber(const int ObjNum){	
	m_numObj=ObjNum;		
}

int Problem::getId() const{ 
	return m_id;     
}

Compare Problem::getOptType(int idx) const{  
	return m_OptMode[idx];   
}

void Problem::resetEvaluations(){            
	m_evals=0;        
}

void Problem::setOptType(const Compare rT,int idx){ 
	if (idx == -1){
		for (auto &i : m_OptMode) i = rT;
	}
	else{
		m_OptMode[idx] = rT;
	}
	      
}

void Problem::setPopInitialMode(PopInitMethod m){
	m_popInitialMode=m;
}
void Problem::setValidateMode(SolutionValidation m){
	m_validationMode=m;
}


bool Problem::isSameType(){ 
	return m_sameType;
}
void Problem::setSameType(bool flag){ 
	m_sameType=flag;
}


void Problem::setProTag(const vector<ProTag> &tag){
	m_tag=tag;
}

void Problem::addProTag(ProTag tag){
	if(!isProTag(tag))	m_tag.push_back(tag);
}
bool Problem::isProTag(ProTag tag){
	auto i=m_tag.begin();
	for(;i!=m_tag.end();++i){
		if(*i==tag) break;
	}
	if(i==m_tag.end()) return false;
	return true;
}

CompareResultFlag Problem::compare(const VirtualEncoding &s1, const VirtualEncoding &s2)const{
	unsigned int better = 0, worse = 0, equal = 0;
	for (unsigned i = 0; i<s1.m_obj.size(); i++){
		if (m_OptMode[i] == MIN_OPT){
			if (s1.m_obj[i]<s2.m_obj[i]) {
				better++;
			}
			else if (s1.m_obj[i]>s2.m_obj[i]) {
				worse++;
			}
			else{
				equal++;
			}
		}
		else{
			if (s1.m_obj[i]>s2.m_obj[i]) {
				better++;
			}
			else if (s1.m_obj[i]<s2.m_obj[i]) {
				worse++;
			}
			else{
				equal++;
			}
		}
	}
	if (better != 0 && better + equal == s1.m_obj.size()) return Compare_Better;
	else if (worse != 0 && worse + equal == s1.m_obj.size()) return Compare_Worse;
	else if (equal == s1.m_obj.size()) return Compare_Equal;
	else return Compare_Non_Dominated;
}
void Problem::copyChanges(const Problem * op, const vector<int> *cd, const vector<int> *co){
	m_evals = op->m_evals;
	m_tevals = op->m_tevals;
	m_cevals = op->m_cevals;

	if (!co){
		for (auto i = 0; i < m_numObj; ++i){
			m_OptMode[i] = op->m_OptMode[(*co)[i]];
		}
	}
}
void Problem::resizeObj(int num){
	m_OptMode.resize(num);
}