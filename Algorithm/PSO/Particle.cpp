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
#include "Particle.h"
#include "../../Problem/ContinuousProblem.h"
#include "../../Global/global.h"

Particle::Particle():Individual(),m_pbest(),\
	m_vel((Global::msp_global->mp_problem.get()!=nullptr)?GET_NUM_DIM:0),\
	m_vMax((Global::msp_global->mp_problem.get()!=nullptr)?GET_NUM_DIM:0){
	initializeVmax();
}
Particle::~Particle(){
	m_vMax.clear();
	m_vMax.clear();
}
void Particle::initializeVmax(){
    for(int i=0;i<GET_NUM_DIM;i++){
	    double l,u;
	    CAST_PROBLEM_CONT->getSearchRange(l,u,i);
		// for static optimization problems
		if(gGetProblemName(Global::ms_curProId).find("FUN_")!=string::npos)
			m_vMax[i].m_max=(u-l)/2;
		else // for dynamic problems
			m_vMax[i].m_max=(u-l)/5;
		m_vMax[i].m_min=-m_vMax[i].m_max;

	}

}

Particle::Particle( const Particle& other):Individual(other),m_pbest(other.m_pbest),m_vel(other.m_vel),m_vMax(other.m_vMax){
    //copy ctor
}

Particle::Particle( const Solution<CodeVReal>& other):Individual(other),	m_pbest(other),m_vel(other.getNumDim()),m_vMax(other.getNumDim()){	
	initializeVmax();
	initializeVelocity();
	m_flag=0;
}
Particle& Particle::operator=(const Particle& rhs){
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    Individual::operator=(rhs);
	m_pbest=rhs.m_pbest;
	m_vel=rhs.m_vel;
	m_vMax=rhs.m_vMax;
    return *this;
}

ReturnFlag Particle::initialize(bool mode){

	ReturnFlag rf= Individual::initialize(mode);
	initializeVelocity();
	m_pbest=self();
	m_flag=0;
	return rf;
}
ReturnFlag Particle::initialize( int ridex, int rid,bool mode){

	ReturnFlag rf= Individual::initialize(ridex,rid,mode);
	initializeVelocity();

	m_pbest=self();
	m_flag=0;
	return rf;
}

ReturnFlag Particle::initialize( int rIdx, int rID,  int rPopsize,bool mode){
	ReturnFlag rf=  Individual::initialize(rIdx,rID,rPopsize,mode);
	initializeVelocity();

	m_pbest=self();
	m_flag=0;
	return rf;
}
void Particle::initialize( const Solution<CodeVReal> &p, int ridex, int rid){

    Individual::initialize(p,ridex,rid);
	initializeVelocity();
	m_pbest=self();
	m_flag=0;

}
ReturnFlag Particle::initialize(const Solution<CodeVReal> &p,double radius,  int ridex, int rid,bool mode){

	ReturnFlag rf= Individual::initialize(p,radius,ridex,rid,mode);
	initializeVelocity();
	m_pbest=self();
	m_flag=0;
	return rf;
}


double Particle::getVelocity(){
	double ve=0;
	for( int i=0;i<GET_NUM_DIM;i++)
		ve+=m_vel[i]*m_vel[i];
	if(ve==0.0) return 0;
	return sqrt(ve);
}
void Particle::initializeVelocity(){
	double u,l;
	for( int i=0;i<GET_NUM_DIM;i++){
	    CAST_PROBLEM_CONT->getSearchRange(l,u,i);
		m_vel[i]=(u-l)*(-0.5+Global::msp_global->mp_uniformAlg->Next())/2;
	}
}

ReturnFlag Particle::move(const Solution<CodeVReal> &lbest, double w, double c1, double c2) {
	double u, l;

	for (int j = 0; j<GET_NUM_DIM; j++) {
		CAST_PROBLEM_CONT->getSearchRange(l, u, j);
		double x = data()[j];
		m_vel[j] = w*m_vel[j] + c1*Global::msp_global->mp_uniformAlg->Next()*(m_pbest.data()[j] - x) + c2*Global::msp_global->mp_uniformAlg->Next()*(lbest.data()[j] - x);

		if (m_vel[j]>m_vMax[j].m_max)	m_vel[j] = m_vMax[j].m_max;
		else if (m_vel[j]<m_vMax[j].m_min)		m_vel[j] = m_vMax[j].m_min;

		data()[j] += m_vel[j];
	}
	self().validate();
	return self().evaluate();
}

ReturnFlag Particle::moveBound( const Solution<CodeVReal> & lbest,double w, double c1, double c2){
	double u,l,cur_x,x;

	for( int j=0;j<GET_NUM_DIM;j++){
		CAST_PROBLEM_CONT->getSearchRange(l,u,j);
		x=(data()[j]);
		if(x==l||x==u){
			double p=Global::msp_global->mp_uniformAlg->Next();
			x=p*(m_pbest.data()[j])+(1-p)*(x);
		}
		cur_x=x;

		double r1=Global::msp_global->mp_uniformAlg->Next();
		double r2=Global::msp_global->mp_uniformAlg->Next();

		m_vel[j]=w*m_vel[j]+c1*r1*((m_pbest.data()[j])-(x))+c2*r2*((lbest.data()[j])-(x));
		if(m_vel[j]>m_vMax[j].m_max)	m_vel[j]=m_vMax[j].m_max;
		else if(m_vel[j]<m_vMax[j].m_min)		m_vel[j]=m_vMax[j].m_min;
		
		x=(x)+m_vel[j];

		while((x)>u||(x)<l){
			x=(cur_x-l<u-cur_x)?cur_x-(cur_x-l)*Global::msp_global->mp_uniformAlg->Next():cur_x+(u-cur_x)*Global::msp_global->mp_uniformAlg->Next();
		} 
		data()[j]=x;
	}
	return self().evaluate();
}


ReturnFlag Particle::NormalMutation(double *avg_v){
	double l,u,x;
	for( int j=0;j<GET_NUM_DIM;j++){
		x=data()[j];
		CAST_PROBLEM_CONT->getSearchRange(l,u,j);
		//if(Global::msp_global->mp_uniformAlg->Next()<1/GET_NUM_DIM)
		x=x+avg_v[j]*Global::msp_global->mp_normalAlg->Next();
		if((x)>u||(x)<l)	x=l+(u-l)*Global::msp_global->mp_uniformAlg->Next();
		data()[j]=x;
	}
	return self().evaluate();
}


void Particle::printToFile(ofstream & out){
        m_pbest.printToFile(out);
}

void Particle::increaseDimension(){

    Individual::increaseDimension();
	m_vel.resize(GET_NUM_DIM);
	m_vMax.resize(GET_NUM_DIM);
    initializeVelocity();
    initializeVmax();
    m_pbest.increaseDimension();

}
void Particle::decreaseDimension(){
    Individual::decreaseDimension();
    m_vel.resize(GET_NUM_DIM);
	m_vMax.resize(GET_NUM_DIM);
    initializeVelocity();
    initializeVmax();
    m_pbest.decreaseDimension();
}
Solution<CodeVReal> & Particle::representative(){
    return m_pbest;
}
const Solution<CodeVReal> & Particle::representative() const{
	return m_pbest;
}
void Particle::updateMemory(){
    Individual::updateMemory();
    m_pbest.evaluate(false);
	if(m_pbest<self()){ 
		Solution<CodeVReal> t=m_pbest;
		m_pbest=self();
		self()=t;
	}
}


void Particle::setVmax(double rMin, double rMax){
	for(int i=0;i<GET_NUM_DIM;i++){
		m_vMax[i].m_max=rMax;
		m_vMax[i].m_min=rMin;

	}

}
void Particle::setVmax(double *rMin, double *rMax){
	for(int i=0;i<GET_NUM_DIM;i++){
	    m_vMax[i].m_max=rMax[i];
		m_vMax[i].m_min=rMin[i];

	}

}
void Particle::printToScreen(){
	Individual::printToScreen();
	cout<<"pbest: ";
	m_pbest.printToScreen();
	cout<<"vel:   ";
	for(auto i=0;i<GET_NUM_DIM;i++){
		cout<<m_vel[i]<<" ";
	}
	cout<<endl;
}

bool Particle::isSame(const Particle &p)const{
	return self().isSame(p.self())&&m_pbest.isSame(p.m_pbest);
}