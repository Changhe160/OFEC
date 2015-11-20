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
// Created: 21 September 2011
// Last modified:
#include "DEIndividual.h"
#include "../../Problem/ContinuousProblem.h"
#include "../../Global/global.h"
DEIndividual::DEIndividual():Individual()
{
    //ctor
}

DEIndividual::~DEIndividual()
{
    //dtor
}

DEIndividual::DEIndividual( const DEIndividual &p):Individual(p){
	m_pv=p.m_pv;
	m_pu=p.m_pu;
}
DEIndividual::DEIndividual( const Solution<CodeVReal> &p):Individual(p){
	m_pv=p;
	m_pu=p;
}
DEIndividual & DEIndividual::operator=(const DEIndividual &rhs){
 if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    Individual::operator=(rhs);
	m_pv=rhs.m_pv;
	m_pu=rhs.m_pu;
	return *this;
}
ReturnFlag DEIndividual::initialize(bool mode){
	ReturnFlag rf= Individual::initialize(-1,-1,mode);
	m_flag=true;
	m_pv=self();
	m_pu=self();
	return rf;
}
ReturnFlag DEIndividual::initialize( int ridex, int rid,bool mode){
	ReturnFlag rf=Individual::initialize(ridex,rid,mode);
	m_flag=true;
    m_pv=self();
	m_pu=self();
	return rf;
}

void DEIndividual::initialize(const Solution<CodeVReal> &p, int idex, int id){

	Individual::initialize(p,idex,id);
	m_flag=true;
    m_pv=self();
	m_pu=self();
}
ReturnFlag DEIndividual::initialize(int rIdx,int rID,int rPopsize,bool mode ){
	ReturnFlag rf=Individual::initialize(rIdx,rID,rPopsize,mode);
	m_flag=true;
	m_pv=self();
	m_pu=self();
	return rf;
}
ReturnFlag DEIndividual::initialize(const Solution<CodeVReal> &p,double radius,  int ridex, int rid,bool mode){
	ReturnFlag rf=Individual::initialize(p,radius,ridex, rid,mode);
	m_flag=true;
    m_pv=self();
	m_pu=self();
	return rf;
}
void DEIndividual::mutate( double F, Solution<CodeVReal> *r1,  Solution<CodeVReal> *r2,  Solution<CodeVReal> *r3, Solution<CodeVReal> *r4,  Solution<CodeVReal> *r5){

	double l,u;
	for(int i=0;i<GET_NUM_DIM;i++){
		CAST_PROBLEM_CONT->getSearchRange(l,u,i);
		m_pv.data()[i]=(r1->data()[i]) +F*((r2->data()[i])-(r3->data()[i]));
		if(r4&&r5) m_pv.data()[i]+=F*((r4->data()[i])-(r5->data()[i]));

		if((m_pv.data()[i])>u){
			m_pv.data()[i]=((r1->data()[i])+u)/2;
		} else if((m_pv.data()[i])<l) {
			m_pv.data()[i]=((r1->data()[i])+l)/2;
		}

	}
}
void DEIndividual::recombine(double CR){
	int I=Global::msp_global->getRandInt(0,GET_NUM_DIM,Program_Algorithm);

	for(int i=0;i<GET_NUM_DIM;i++){
		double p=Global::msp_global->mp_uniformAlg->Next();
		if(p<=CR||i==I)     m_pu.data()[i]=m_pv.data()[i];
		else m_pu.data()[i]=data()[i];
	}

}
ReturnFlag DEIndividual::select(){

    ReturnFlag rf=m_pu.evaluate();
    if(m_pu>self()) self()=m_pu;
	return rf;
}
void DEIndividual::increaseDimension(){
    Individual::increaseDimension();
	m_pv.increaseDimension();
	m_pu.increaseDimension();
}
void DEIndividual::decreaseDimension(){
	Individual::decreaseDimension();
    m_pv.decreaseDimension();
	m_pu.decreaseDimension();
}
void DEIndividual::printToScreen(){
	Individual::printToScreen();
	cout<<"pv  : ";
	m_pv.printToScreen();
	cout<<"pu  : ";
	m_pu.printToScreen();
}

