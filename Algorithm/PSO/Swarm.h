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
// Created: 21 September 2011
// Last modified: 12 Dec. 2014

#ifndef SWARM_H
#define SWARM_H
#include "../PopulationCont.h"

template <typename ED,typename TypeParticle>
class Swarm: public PopulationCont<ED,TypeParticle>{
protected:
	double m_C1,m_C2;                       //accelerators
	double m_W,m_maxW, m_minW;                // inertia weight
	//May 10, 2016
	vector<vector<bool>> m_link;
public:
	virtual ~Swarm(void){}
	virtual void initialize( bool mode=true);

	explicit Swarm();
	Swarm(int rPopsize,bool mode);
	Swarm(const Swarm<ED,TypeParticle> &s);
	Swarm(Group<ED,TypeParticle> &g);
	Swarm(const TypeParticle & center, double radius, int rPopsize,bool mode);

	Swarm &operator= (const Swarm &s);           
	void initializePara(double rw,double rc1,double rc2,bool write2File=true,double rmaxw=0,double rminw=0);
	void setW(const double rw);
	double getW() const { return m_W; }
	void setC1(const double rc1);
	double getC1() const { return m_C1; }
	void setC2(const double rc2);
	double getC2() const { return m_C2; }

	virtual void reInitialize(bool clearOldBest,bool mode=true);
	void reInitialize(const TypeParticle & center,bool mode=true);
	double getAvgVelocity(); 
protected:
	//May 10, 2016
	virtual void setNeibourhood(); //gbest as default 
	virtual Solution<ED>& neighborBest(int);
};


template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle>::Swarm():PopulationCont<ED,TypeParticle>(){
	// for static optimization problems
	if(gGetProblemName(Global::ms_curProId).find("FUN_")!=string::npos)
		initializePara(0.7298f,1.496f,1.496f,true);
	else // for dynamic problems
		initializePara(0.6f,1.7f,1.7f,false,0.6f,0.3f);
}

template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle>::Swarm(const int rPopsize, bool mode):PopulationCont<ED,TypeParticle>(rPopsize,mode){
	// for static optimization problems
    initialize();

}
template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle>::Swarm(const Swarm &s):PopulationCont<ED,TypeParticle>(s){
	m_W=s.m_W;
	m_C1=s.m_C1;
	m_C2=s.m_C2;
	m_maxW=s.m_maxW;
	m_minW=s.m_minW;
	m_link = s.m_link;
}
template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle>::Swarm(Group<ED,TypeParticle> &g):PopulationCont<ED,TypeParticle>(g){

	initialize();
}
template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle>::Swarm(const TypeParticle & center, double radius,const int rPopsize,bool mode):PopulationCont<ED,TypeParticle>(center,radius,rPopsize,mode){
	initialize();
}
template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::initialize(bool mode){

    // for static optimization problems
	if(gGetProblemName(Global::ms_curProId).find("FUN_")!=string::npos)
		initializePara(0.7298f,1.496f,1.496f);
	else // for dynamic problems
		initializePara(0.6f,1.7f,1.7f,0,0.6f,0.3f);

	//**************** July 20 2012*****************//
	// constrain particles' max velocity wth the swarm radius
	if(Global::ms_curProId==Global::msm_alg["ALG_AMSO"])
		for(int i=0;i<this->m_popsize;i++){
			this->m_pop[i]->setVmax(-this->m_initialRadius,this->m_initialRadius);
		}
	//**************** July 20 2012*****************//
	m_link.resize(this->m_popsize);
	for (auto &i : m_link) i.resize(this->m_popsize);

}
template <typename ED,typename TypeParticle>
Swarm<ED,TypeParticle> &Swarm<ED,TypeParticle>::operator= (const Swarm &s){
	if(this==&s) return *this;
	PopulationCont<ED,TypeParticle>::operator=(s);

	m_W=s.m_W;
	m_C1=s.m_C1;
	m_C2=s.m_C2;
	m_maxW=s.m_maxW;
	m_minW=s.m_minW;
	m_link = s.m_link;
	return *this;
}

template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::initializePara(double rw,double rc1,double rc2,bool write2File,double rmaxw,double rminw){
    if(write2File){
	setW(rw);
	setC2(rc2);
	setC1(rc1);
    }else{
        m_W=rw; m_C1=rc1;m_C2=rc2;
    }
	m_maxW=rmaxw;
	m_minW=rminw;
}

template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::reInitialize(bool clearOldBest,bool mode){
	for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->initialize(i,i+1,mode); // index, ID
	if(this->ifTerminating()) return;
	
	this->computeCenter();
	this->computeInitialRadius();
	this->findBest();
	this->findWorst();
	if(clearOldBest) this->m_best.clear();

	for(unsigned i=0;i<this->m_bestIdx.size();i++){
		this->updateBestArchive(this->m_pop[this->m_bestIdx[i]]->representative());
	} 
}

template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::reInitialize(const TypeParticle & center, bool mode){
	PopulationCont<ED,TypeParticle>::initialize(center,this->m_initialRadius,mode);
}


template <typename ED,typename TypeParticle>
double Swarm<ED,TypeParticle>::getAvgVelocity(){
	double avg=0;
	if(this->m_popsize<1) return 0;
	for(int i=0;i<this->m_popsize;i++) avg+=this->m_pop[i]->getVelocity();

	return avg/this->m_popsize;

}

template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::setW(const double rw){
	m_W=rw;

	size_t start, end;
    start=this->m_algPar.str().find("Inertia weight:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;

    ss<<"Inertia weight:"<<m_W<<"; ";
	if(start!=string::npos){
		string result=this->m_algPar.str();
		result.replace(start,end-start+1, ss.str());
		 this->m_algPar.str(result);
	}else this->m_algPar<<ss.str();
}
template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::setC1(const double rc1){
	m_C1=rc1;

	size_t start, end;
    start=this->m_algPar.str().find("Accelerate constant1:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Accelerate constant1:"<<m_C1<<"; ";
	if(start!=string::npos){
		string result=this->m_algPar.str();
		result.replace(start,end-start+1, ss.str());
		 this->m_algPar.str(result);
	}else this->m_algPar<<ss.str();


}
template <typename ED,typename TypeParticle>
void Swarm<ED,TypeParticle>::setC2(const double rc2){
	m_C2=rc2;

	size_t start, end;
    start=this->m_algPar.str().find("Accelerate constant2:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Accelerate constant2:"<<m_C2<<"; ";
	if(start!=string::npos){
		string result=this->m_algPar.str();
		result.replace(start,end-start+1, ss.str());
		 this->m_algPar.str(result);
	}else this->m_algPar<<ss.str();

}

template <typename ED, typename TypeParticle>
void Swarm<ED, TypeParticle>::setNeibourhood() {
	for (auto &i : m_link) {
		for (int j = 0; j < i.size(); ++j) {
			i[j] = true;
		}
	}
}
template <typename ED, typename TypeParticle>
Solution<ED>& Swarm<ED, TypeParticle>::neighborBest(int idx) {
	return *this->m_best[0];
}
#endif // SWARM_H
