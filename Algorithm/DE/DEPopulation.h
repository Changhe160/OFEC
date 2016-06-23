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

/*Storn, R. and Price, K. (1997), "Differential Evolution - A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces",
Journal of Global Optimization, 11, pp. 341-359*/
#ifndef DEPOPULATION_H
#define DEPOPULATION_H
#include "../PopulationCont.h"

template<typename ED,typename TypeDEIndi>
class DEPopulation: public PopulationCont<ED,TypeDEIndi>{
protected:
    double  m_F, m_CR;
    DEMutationStratgy m_mutStrategy;
public:
    DEPopulation();
    virtual ~DEPopulation(){}

    DEPopulation(int rPopsize,bool mode=true);
	DEPopulation(const DEPopulation<ED,TypeDEIndi> &s);
	DEPopulation(Group<ED,TypeDEIndi> &g);
	DEPopulation(const Solution<CodeVReal> & center, double radius, int rPopsize,bool mode=true);
	void setMutationStrategy(DEMutationStratgy rS);
    DEPopulation & operator=(const DEPopulation & s);
    virtual void mutate(const int idx);
    void setParmeter(double cr, double f);
	void defaultParameter();
	virtual void reInitialize(bool clearOldBest,bool mode);	
protected:
	virtual void select(int base, int number, vector<int>&);
	ReturnFlag evolve();
};
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi>::DEPopulation():PopulationCont<ED,TypeDEIndi>(),m_F(0.5),m_CR(0.1),m_mutStrategy(DE_rand_1){
		defaultParameter();
}
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi>::DEPopulation(int rPopsize,bool mode):PopulationCont<ED,TypeDEIndi>(rPopsize,mode){
	defaultParameter();
}
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi>::DEPopulation(const DEPopulation<ED,TypeDEIndi> &s):PopulationCont<ED,TypeDEIndi>(s),m_F(s.m_F ),m_CR(s.m_CR),m_mutStrategy(s.m_mutStrategy){

}
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi>::DEPopulation(Group<ED,TypeDEIndi> &g):PopulationCont<ED,TypeDEIndi>(g){
	defaultParameter();
}
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi>::DEPopulation(const Solution<CodeVReal> & center, double radius, int rPopsize,bool mode):PopulationCont<ED,TypeDEIndi>(center,radius,rPopsize,mode){
	defaultParameter();
}
template<typename ED,typename TypeDEIndi>
void DEPopulation<ED,TypeDEIndi>::setMutationStrategy(DEMutationStratgy rS){
    m_mutStrategy=rS;
}
template<typename ED,typename TypeDEIndi>
DEPopulation<ED,TypeDEIndi> & DEPopulation<ED,TypeDEIndi>::operator=(const DEPopulation & s){
    if(this==&s) return *this;
	DEPopulation<ED,TypeDEIndi>::operator=(s);
    m_CR=s.m_CR;
    m_F=s.m_F;
    m_mutStrategy=s.m_mutStrategy;
    return *this;
}
template<typename ED,typename TypeDEIndi>
void DEPopulation<ED,TypeDEIndi>::mutate(const int idx){
	vector<int> ridx;
    switch(m_mutStrategy){
        case DE_rand_1:
			select(idx, 3, ridx);
            this->m_pop[idx]->mutate(m_F,&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self(),&this->m_pop[ridx[2]]->self());
            break;
        case DE_best_1:
			select(idx, 2, ridx);
            this->m_pop[idx]->mutate(m_F,&this->m_best[0]->self(),&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self());
            break;
        case DE_targetToBest_1:
			select(idx, 2, ridx);
            this->m_pop[idx]->mutate(m_F,&this->m_pop[idx]->self(),&this->m_best[0]->self(),&this->m_pop[idx]->self(),&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self());
            break;
        case DE_best_2:
			select(idx, 4, ridx);
            this->m_pop[idx]->mutate(m_F,&this->m_best[0]->self(),&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self(),&this->m_pop[ridx[2]]->self(),&this->m_pop[ridx[3]]->self());
            break;
        case DE_rand_2:
			select(idx, 5, ridx);
            this->m_pop[idx]->mutate(m_F,&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self(),&this->m_pop[ridx[2]]->self(),&this->m_pop[ridx[3]]->self(),&this->m_pop[ridx[4]]->self());
            break;
		case DE_randToBest_1:
			select(idx, 3, ridx);
			this->m_pop[idx]->mutate(m_F,&this->m_pop[ridx[0]]->self(),&this->m_best[0]->self(),&this->m_pop[ridx[0]]->self(),&this->m_pop[ridx[1]]->self(),&this->m_pop[ridx[2]]->self());
			break;
		case DE_targetToRand_1:
			select(idx, 3, ridx);
			this->m_pop[idx]->mutate(m_F,&this->m_pop[idx]->self(),&this->m_pop[ridx[0]]->self(),&this->m_pop[idx]->self(),&this->m_pop[ridx[1]]->self(),&this->m_pop[ridx[2]]->self());
			break;
    }
}
template<typename ED,typename TypeDEIndi>
ReturnFlag DEPopulation<ED,TypeDEIndi>::evolve(){
    if(this->m_popsize<5){
		throw myException("the population size cannot be smaller than 5@DEPopulation<ED,TypeDEIndi>::evolve()");
    }
	ReturnFlag r_flag=Return_Normal;
    for(int i=0;i<this->m_popsize;i++){
        mutate(i);
        this->m_pop[i]->recombine(m_CR);
    }

    this->updateIDnIndex();
    for(int i=0;i<this->m_popsize;i++){
        r_flag=this->m_pop[i]->select();
		this->updateBestArchive(this->m_pop[i]->self());
		if(r_flag!=Return_Normal) break;

    }
	if(r_flag==Return_Normal) {
		this->m_center=*this->m_best[0];
		this->updateCurRadius();
		this->m_iter++;
	}
    return r_flag;
}
template<typename ED,typename TypeDEIndi>
void DEPopulation<ED,TypeDEIndi>::setParmeter(const double cr, const double f){
    m_CR=cr;
    m_F=f;
}
template<typename ED,typename TypeDEIndi>
void DEPopulation<ED,TypeDEIndi>::defaultParameter(){
	m_CR=0.6;
	m_F=0.5;
	m_mutStrategy=DE_best_1;
}
template<typename ED,typename TypeDEIndi>
void DEPopulation<ED,TypeDEIndi>::reInitialize(bool clearOldBest,bool mode){
	for(int i=0;i<this->m_popsize;i++){
		this->m_pop[i]->initialize(mode);
		this->m_pop[i]->m_flag=true;
	}
	if(this->ifTerminating()) return;
	this->findBest();
	if(clearOldBest) this->m_best.clear();
	for(unsigned i=0;i<this->m_bestIdx.size();i++){
			this->updateBestArchive(this->m_pop[this->m_bestIdx[i]]->representative());
	}
}
template<typename ED, typename TypeDEIndi>
void DEPopulation<ED, TypeDEIndi>::select(int base, int number, vector<int>& result){
	vector<int> candidate;
	for (int i = 0; i < this->m_popsize; i++) {
		if (base != i) candidate.push_back(i);
	}
	for (int i = 0; i < number; ++i) {
		int idx=Global::msp_global->getRandInt(0, candidate.size());
		result.push_back(candidate[idx]);
		candidate.erase(candidate.begin() + idx);
	}
}
#endif // DEPOPULATION_H


