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

#ifndef MULTI_POPULATION_H
#define MULTI_POPULATION_H
#include "Population.h"
template<typename TypePop>
class MultiPopulation{
protected:
	 vector<unique_ptr<TypePop> > m_subPop;   // populations, warning: it currently can only contain populations belong to the same type
	 int m_maxSubSize;			//max size for each population
	 int m_convergedPops;
public:  
	virtual ~MultiPopulation();
	MultiPopulation():m_maxSubSize(0),m_convergedPops(0){}
	MultiPopulation(int size):m_maxSubSize(size),m_convergedPops(0){}
	MultiPopulation(const int num,const int subSize,bool mode, const int runId);

    virtual void printPops(ofstream &out,Problem *pro=0);
  
	void updateMemoryAll();
	int findWorstPop();
	int findBestPop(int pos=0,int flagIdx=-1,bool mode=true);
	void increaseDimensionAll();
	void decreaseDimensionAll();
	void changeDimensionAll();

	void ROOTBestAll();
	void addPopulation(TypePop & w);
	void deletePopulation(unsigned index);
	unsigned getNumPops();
	const unique_ptr<TypePop> & operator[](const int i) const;
	unique_ptr<TypePop> & operator[](const int i);
	virtual void handleReturnFlagAll(ReturnFlag f);
	int computePeaksFoundAll();
	void setMaxSubSize(unsigned size);
protected:
	void setPopID(TypePop & w);
};

template<typename TypePop>
	MultiPopulation<TypePop>::~MultiPopulation(){
	for(auto &pop:m_subPop){
		pop.reset();
	}
	m_subPop.clear();
}

template<typename TypePop>
MultiPopulation<TypePop>::MultiPopulation(const int num,const int subSize,bool mode, const int runId):m_maxSubSize(subSize),m_convergedPops(0){
	for(int i=0;i<num;i++)		m_subPop.push_back(new TypePop(subSize,mode,runId));
}

template<typename TypePop>
void MultiPopulation<TypePop>::printPops(ofstream &out,Problem *pro){	
	for(auto &pop:m_subPop){
			pop->printPopToFile(out);	
	}
}

template<typename TypePop>
void MultiPopulation<TypePop>::updateMemoryAll(){

	for(auto &pop:m_subPop){
			pop->updateMemory();	
	}

}

template<typename TypePop>
void MultiPopulation<TypePop>::setPopID(TypePop & w){
	int maxId=-1;
	for(auto &i:m_subPop){
		if(i->getID()>maxId) maxId=i->getID();
	}
	if(maxId==-1) w.setID(1);
	else{
		vector<int> vid(maxId);
		for(auto i=vid.begin();i!=vid.end();++i){
			*i=i-vid.begin()+1;
		}
		for(auto &i:m_subPop){
			decltype (vid.size()) j=0;
			for(;j<vid.size();j++){
				if(i->getID()==vid[j]) break;
			}
			if(j<vid.size()) vid.erase(vid.begin()+j);
		}
		if(vid.empty()) w.setID(maxId+1);
		else w.setID(vid[0]);
	}
}

template<typename TypePop>
int MultiPopulation<TypePop>::findWorstPop(){
	//return the pop index with the worst best individual 
	if(m_subPop.empty()) return -1;
    
	typename vector<unique_ptr<TypePop> >::iterator it=m_subPop.begin();
	int idx=0,i=0;
	for(i=0;it!=m_subPop.end();++it,++i) {
		if((*it)->m_popsize==0) continue;
		if((*it)<(*(it+idx))){
			idx=i;
		}
	}

	return idx;
}

template<typename TypePop>
int MultiPopulation<TypePop>::findBestPop(int pos,int flagIdx,bool mode){
	//return the pop index with the best best individual with m_flag[flagIdx]=mode between m_subPop.begin()+pos, and m_subPop.end() 
	if(m_subPop.empty()) return -1;
	typename vector<unique_ptr<TypePop> >::iterator it=m_subPop.begin()+pos;
	typename vector<unique_ptr<TypePop> >::iterator bestIt=it;
	for(;it!=m_subPop.end();++it) {
		if((*it)->m_popsize==0) continue;
		if(flagIdx!=-1&&(*it)->m_flag[flagIdx]!=mode) continue;
		if(*(*it)>*(*(bestIt))){
			bestIt=it;
		}
	}
	return bestIt-m_subPop.begin();
}

template<typename TypePop>
void MultiPopulation<TypePop>::increaseDimensionAll(){
	for(auto &pop:m_subPop){
		pop->increaseDimension();
	}
}
template<typename TypePop>
void MultiPopulation<TypePop>::decreaseDimensionAll(){
	for(auto &pop:m_subPop){
		pop->decreaseDimension();
	}
}

template<typename TypePop>
void MultiPopulation<TypePop>::changeDimensionAll(){
	for(auto &pop:m_subPop){
		pop->changeDimension();
	}
}

template<typename TypePop>
void MultiPopulation<TypePop>::ROOTBestAll(){
	if(mROOT::getROOT()){
		for(unsigned int j=0;j<m_subPop.size();j++) for(auto& j:m_subPop[j]->m_best)mROOT::getROOT()->record(j);
	}
}

template<typename TypePop>
void MultiPopulation<TypePop>::addPopulation(TypePop & w){
	setPopID(w);
	m_subPop.push_back(move(unique_ptr<TypePop>(&w)));	
}

template<typename TypePop>
void MultiPopulation<TypePop>::deletePopulation(unsigned index){
	if(index<0||index>=m_subPop.size())  throw myException("no subswarm of index@deletePopulation(unsigned index)");
	//m_subPop[index].reset();
	m_subPop.erase(m_subPop.begin()+index);
}

template<typename TypePop>
unsigned MultiPopulation<TypePop>::getNumPops(){
	return m_subPop.size();
}

template<typename TypePop>
const unique_ptr<TypePop> & MultiPopulation<TypePop>::operator[](const int i) const{
	return m_subPop[i];
}

template<typename TypePop>
unique_ptr<TypePop> & MultiPopulation<TypePop>::operator[](const int i){
	return m_subPop[i];
}
	
template<typename TypePop>
void MultiPopulation<TypePop>::handleReturnFlagAll(ReturnFlag f){
	for(unsigned i=0;i<m_subPop.size();i++){
		m_subPop[i]->handleReturnFlag(f);
	};
}


template<typename TypePop>
void MultiPopulation<TypePop>::setMaxSubSize(unsigned size){
	m_maxSubSize=size;
}
#endif