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
// Last modified: 12 Dec. 2014
#ifndef POPULATION_H
#define POPULATION_H
#include "Algorithm.h"
#include "../Utility/Clustering/Group.h"
#include "../Problem/DOP/DynamicContinuous.h"
#include "../Measure/mROOT.h"
#include "../Measure/mMultiModal.h"
#include "../Measure/mSingleObj.h"
#include "../Problem/optimum.h"
#include "../Problem/FunctionOpt/BenchmarkFunction.h"
template<typename>  class MultiPopulation;

template <typename ED, typename TypeIndi>
class Population: public  Algorithm{
	template<typename> friend class MultiPopulation;
protected:
	int m_popsize;                              // swarm size
	int m_maxID;                                // the maximum ID of the population
	vector<unique_ptr<TypeIndi>> m_pop;			// the population, the number of indis can changes
	int m_evoNum;                               // the number of generations
	int m_popID;                                // the ID of the populaton
	vector<int> m_bestIdx,m_worstIdx;                   // the indice of the best and worst individual(warning: for non-dominate population, they are equal to -1)
	vector<bool> m_flag;                                // for use in some cases, e.g., whether it found an optimum solution or not, ....
	vector<unique_ptr<Solution<ED>>> m_best;              // the best one in the population, notice: can be non-donimated solutions for multi-objective optimization problems	
	vector<int> m_orderList;							// indices of ordered individuals
	int m_curRank=0;
protected:
	void copy_(const Population<ED, TypeIndi> &s);
	virtual ReturnFlag initialize(bool isInitialized,bool mode, bool clearOldBest);
	void update(const int num);
	inline void updateIndex();
	inline void updateIDnIndex();
	inline void updateMaxID();
	bool updateBestArchive(const Solution<ED> & chr);	
	const Solution<ED> & getNearestBest(const Solution<ED> &chr);
public:
	virtual ~Population(void);
	explicit Population():Algorithm(-1,string()),m_popsize(0),m_maxID(0),m_evoNum(0),m_popID(0),m_orderList(0),m_best(0),m_pop(0){}	
	Population(const int rPopsize);
	Population(const int rPopsize,bool mode);
	Population(const Population &s);
	//transfer ownership to m_pop
	Population(Group<ED,TypeIndi> &g);
	Population & operator= (const Population &s);
	 const vector<int>& findBest(void);
	const vector<int>& findWorst(void);
	virtual void add( Group<ED,TypeIndi> &g);
	virtual void add(TypeIndi *p,bool tranship=true);
	virtual void add( Population<ED, TypeIndi> &s);
	virtual ReturnFlag add( int num,  bool mode);
  
	virtual void add(vector<TypeIndi*> &indis);
	virtual void add( vector<unique_ptr<TypeIndi>> &indis);
	virtual void remove( int num,const int *id=0);
    virtual void remove( const  vector<int> &id);            
	virtual void increaseDimension();
	virtual void decreaseDimension();
	void changeDimension();
	void printPopToFile(ofstream &out);

	virtual void updateMemory();
	void sort(bool mode, bool representative);
	void rank(); 
	inline double getMaxObj(const int idx);
	inline double getMinObj(const int idx);
	void checkOverCrowd(int subSize);

	void ROOTBest(vector<int> idx);
	void ROOTBest();
	const vector<unique_ptr<TypeIndi>> & getPop();
	inline const int getPopSize()const;
	void handleReturnFlag(ReturnFlag f);
	virtual bool isConverged();
	TypeIndi* operator[](const int i);
	const TypeIndi* operator[](const int i)const;
	bool operator >(const Population<ED, TypeIndi>& p);
	bool operator <(const Population<ED, TypeIndi>& p);
	inline const vector<unique_ptr<Solution<ED>>>& getBest()const;
	void printBest2Screen();
	inline int getID();
	inline void setID(int id);
	inline vector<bool> getFlag() const;
	inline int getEvoNum() const;
	inline int getCurRank();
	void resize(int size);
	double rank(const Solution<ED> & i, bool mode = true){
		return rank(i.data(), mode);
	}
	double rank(const ED & i, bool mode = true);
	void clearBestArchive();
	#ifdef OFEC_DEMON
	void startRankThread(vector<int> &rank_, vector<int> &count, vector<vector<int> >& cset);
	static int rankThread(vector<int> &tsk, const Population<ED, TypeIndi> &p, vector<int> &rank_, vector<int> &count, vector<vector<int> >& cset);
	#endif
};

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::copy_(const Population<ED, TypeIndi> &s){

	for(int i=0;i<m_popsize;i++){
		*m_pop[i]=*s.m_pop[i];
	}
	if(m_best.size()<s.m_best.size()){
		for(unsigned i=m_best.size();i<s.m_best.size();i++){
			m_best.push_back(move(unique_ptr<Solution<ED>>( new Solution<ED>(*s.m_best[i]) )));
		}
	}else if(m_best.size()>s.m_best.size()){
		m_best.erase(m_best.begin()+s.m_best.size(),m_best.end());
	}
	for(unsigned i=0;i<m_best.size();i++){
		*m_best[i]=*s.m_best[i];
	}
		
	m_maxID=s.m_maxID;
	m_evoNum=s.m_evoNum;
	m_popID=s.m_popID;
	m_flag=s.m_flag;
		
	m_orderList=s.m_orderList;
	m_bestIdx=s.m_bestIdx;
	m_worstIdx=s.m_worstIdx;

	m_curRank = s.m_curRank;
}

template <typename ED, typename TypeIndi>
ReturnFlag Population<ED, TypeIndi>::initialize(bool isInitialized,bool mode, bool clearOldBest){
	ReturnFlag rf=Return_Normal;
	if(!isInitialized){	
		for(int i=0;i<m_popsize;i++){
			rf=m_pop[i]->initialize(i,i+1,m_popsize,mode); // index, ID
			if(rf==Return_Terminate) return rf; 
			else if(rf==Return_Change){
				stringstream ss;
				ss<<"an environmental change occurs during population initialization at evals "<<Global::msp_global->mp_problem->getEvaluations()<<" mode "<<mode<<" clearOldBest "<<clearOldBest<<" @ Population<ED, TypeIndi>::initialize(bool ,bool , bool)";
				throw myException(ss.str().c_str());
			} 
		}
	}
	m_evoNum=0;
	findWorst();
	findBest();
	if(clearOldBest) m_best.clear();
	for(unsigned i=0;i<m_bestIdx.size();i++){
		updateBestArchive(m_pop[m_bestIdx[i]]->representative());
	}
	return rf;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::update(const int num){
	for(int i=m_popsize-num;i<m_popsize;i++){
		m_maxID++;
		m_pop[i]->m_index=i; 
		m_pop[i]->m_id=m_maxID;
	}
		
	findWorst();
	findBest();
	for(auto &i:m_bestIdx){
		updateBestArchive(m_pop[i]->representative());
	}
	Global::msp_global->m_totalNumIndis+=num;
	m_orderList.resize(m_popsize);
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::updateIndex(){
	for(int i=0;i<m_popsize;i++)  m_pop[i]->m_index=i;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::updateIDnIndex(){
		for(int i=0;i<m_popsize;i++){
		m_pop[i]->m_index=i;
		m_pop[i]->m_id=i+1;
	}

}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::updateMaxID(){
	int max=0;
	for(int i=0;i<m_popsize;i++)
		if(max<m_pop[i]->m_id) max=m_pop[i]->m_id;
	m_maxID=max;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::printBest2Screen(){
	for(auto i=m_best.begin();i!=m_best.end();i++){
		(*i)->printToScreen();
	}
}

template <typename ED, typename TypeIndi>
bool Population<ED, TypeIndi>::updateBestArchive(const Solution<ED> & chr){
	bool first=true;
	// check dominated case
	for(auto i=m_best.begin();i!=m_best.end();i++){
		if(first && **i<chr){
			**i=chr;
			first=false;
		}else if(!first && **i<chr){
			i=m_best.erase(i);
			i--;
		}
	}
	if(!first) return false;
	//check equal case	
	for(auto i=m_best.begin();i!=m_best.end();i++){
		if(**i==chr&&!((*i)->isSame(chr))){
			m_best.push_back(move(unique_ptr<Solution<ED>>(new Solution<ED>(chr))));
		//	printBest2Screen();
			return true;
		}
	}
	//check non-dominated case	
	for(auto i=m_best.begin();i!=m_best.end();i++){
		if(!(**i!=chr)) return false;
	}
	m_best.push_back(move(unique_ptr<Solution<ED>>(new Solution<ED>(chr))));
	return true;	
}

template <typename ED, typename TypeIndi>
const Solution<ED> & Population<ED, TypeIndi>::getNearestBest(const Solution<ED> &chr){
	if(m_best.size()==0) throw myException("size of m_best=0 @getNearestBest(const Solution &chr)");
	if(m_best.size()==1) return *m_best[0];

	double dis=chr.getDistance(*m_best[0]);
	unsigned idx=0;
	for(unsigned i=1;i<m_best.size();i++){
		double dis_=chr.getDistance(*m_best[i]);
		if(dis>dis_){
			dis=dis_;
			idx=i;
		}	
	}
	return *m_best[idx];
}

template <typename ED, typename TypeIndi>
Population<ED, TypeIndi>::~Population(void){
	bool check[Loki::SuperSubclass<Solution<ED>,TypeIndi>::value]; // TypeIndi must be derived from Solution
	m_pop.clear();
	m_bestIdx.clear();
	m_worstIdx.clear();                   
	m_flag.clear();                                
	m_best.clear();             
	m_orderList.clear();
	if(Global::msp_global.get())	Global::msp_global->m_totalNumIndis-=m_popsize;	
}

template <typename ED, typename TypeIndi>
Population<ED, TypeIndi>::Population(const int rPopsize):Algorithm(-1,string()),m_popsize(rPopsize),m_maxID(rPopsize),m_evoNum(0),m_popID(0),\
	m_flag(1,false),m_best(),m_pop(rPopsize),m_orderList(rPopsize){
	for(auto &i:m_pop) i=move(unique_ptr<TypeIndi>(new TypeIndi()));
	Global::msp_global->m_totalNumIndis+=m_popsize;
}

template <typename ED, typename TypeIndi>
Population<ED, TypeIndi>::Population(const int rPopsize,bool mode):Algorithm(-1,string()),m_popsize(rPopsize),m_maxID(rPopsize),m_evoNum(0),m_popID(0),\
	m_flag(1,false),m_best(),m_pop(rPopsize),m_orderList(rPopsize){
	for(auto &i:m_pop) i=move(unique_ptr<TypeIndi>(new TypeIndi()));
	initialize(false,mode,false);
	Global::msp_global->m_totalNumIndis+=m_popsize;
}

template <typename ED, typename TypeIndi>
Population<ED, TypeIndi>::Population(const Population<ED, TypeIndi> &s):Algorithm(-1,string()),m_popsize(s.m_popsize),m_pop(s.m_popsize),m_orderList(s.m_orderList),m_best(s.m_best.size()){	
		
	for(int i=0;i<m_popsize;i++){
		m_pop[i]=move(unique_ptr<TypeIndi>(new TypeIndi(*s.m_pop[i]))) ;
	}
	for(unsigned i=0;i<m_best.size();i++){
		m_best[i]=move(unique_ptr<Solution<ED>>(new Solution<ED>(*s.m_best[i])));
	}

	m_maxID=s.m_maxID;
	m_evoNum=s.m_evoNum;
	m_popID=s.m_popID;
	m_flag=s.m_flag;
	m_bestIdx=s.m_bestIdx;
	m_worstIdx=s.m_worstIdx;
	m_curRank = s.m_curRank;
	Global::msp_global->m_totalNumIndis+=m_popsize;
}

//transfer ownership to m_pop
template <typename ED, typename TypeIndi> 
Population<ED, TypeIndi>::Population( Group<ED,TypeIndi> &g):Algorithm(-1,string()),m_popsize(g.getSize()),m_maxID(g.getSize()),m_evoNum(0),m_popID(0),m_flag(1,false),\
	m_best(g.getBest().size()),m_pop(g.getSize()),m_orderList(g.getSize()){
		
	for(int i=0;i<m_popsize;i++){
		m_pop[i]=move(unique_ptr<TypeIndi>(g[i]));
		m_pop[i]->m_index=i;
		m_pop[i]->m_id=i+1;
	}
	for(unsigned i=0;i<m_best.size();i++){
		m_best[i]=move(unique_ptr<Solution<ED>>(new Solution<ED>(g.getBest()[i])));
	}

	m_evoNum=0;
	findWorst();
	Global::msp_global->m_totalNumIndis+=m_popsize;	
}

template <typename ED, typename TypeIndi>
Population<ED, TypeIndi> & Population<ED, TypeIndi>::operator= (const Population<ED, TypeIndi> &s){
	if(this==&s) return *this;
	if(m_popsize!=s.m_popsize){	
		throw myException("the size of two populations must be the same in assignment operator@Population & operator= (const Population &s)"); return *this;
	}
	Algorithm::operator=(s);

	copy_(s);
	return *this;
}

template <typename ED, typename TypeIndi>
const vector<int>& Population<ED, TypeIndi>::findBest(void){
	m_bestIdx.clear();
	if(m_popsize<1){
		return m_bestIdx;
	}
	vector<bool> flag(m_popsize,true);

	for(int j=0;j<m_popsize;j++){
		for(int i=0;i<m_popsize;i++){
			if(i==j||!flag[j]||!flag[i]) continue;
			if(*m_pop[j]>*m_pop[i]){
				flag[i]=false;
			}
		}
	}
	for(int i=0;i<m_popsize;i++){ if(flag[i]) m_bestIdx.push_back(i);}
	return m_bestIdx;
}

template <typename ED, typename TypeIndi>
const vector<int>& Population<ED, TypeIndi>::findWorst(void){	
	m_worstIdx.clear();
	if(m_popsize<1){
		return m_worstIdx;
	}
	vector<bool> flag(m_popsize,true);
	for(int j=0;j<m_popsize;j++){
		for(int i=0;i<m_popsize;i++){
			if(i==j||!flag[j]||!flag[i]) continue;
			if(*m_pop[j]<*m_pop[i]){
				flag[i]=false;
			}
		}
	}
	for(int i=0;i<m_popsize;i++){ if(flag[i]) m_worstIdx.push_back(i);}
		
	return m_worstIdx;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::add(Group<ED,TypeIndi> &g){
	//transfer ownership from g.m_member to m_pop
	m_popsize+=g.getSize();
	for(int i=0;i<g.getSize();i++){
		m_pop.push_back(move(unique_ptr<TypeIndi>(g[i])));
	}
	update(g.getSize());
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::add(TypeIndi *p,bool tranship){
	++m_popsize;
	if(tranship)	m_pop.push_back(move(unique_ptr<TypeIndi>(p)));
	else m_pop.push_back(move(unique_ptr<TypeIndi>(new TypeIndi(*p))));
	p=0;
	update(1);
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::add( Population<ED, TypeIndi> &s){
	//transfer ownership from s.m_pop to m_pop
	m_popsize+=s.m_popsize;
	for(auto& i:s.m_pop){
		m_pop.push_back(move(i));
	}
	update(s.m_popsize);
		
}

template <typename ED, typename TypeIndi>
ReturnFlag Population<ED, TypeIndi>::add( int num,  bool mode){
	ReturnFlag rf=Return_Normal;
	//m_popsize+=num;
	int count=0;
	for(int i=0;i<num;i++){
		m_pop.push_back(move(unique_ptr<TypeIndi>(new TypeIndi())));
		rf=m_pop[m_popsize]->initialize(m_popsize,m_popsize+1,m_popsize+1,mode);
		++m_popsize;
		++count;
		if(rf!=Return_Normal) break;
	}
	update(count);
	return rf;
}
  
template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::add(vector<TypeIndi*> &indis){
	//transfer owership here 
	m_popsize+=indis.size();
	for(unsigned i=0;i<indis.size();i++){
		m_pop.push_back(move(unique_ptr<TypeIndi>(indis[i])));
		indis[i]=0;
	}
	update(indis.size());
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::add( vector<unique_ptr<TypeIndi>> &indis){
	//transfer owership here 
	m_popsize+=indis.size();
	for(unsigned i=0;i<indis.size();i++){
		m_pop.push_back(move(indis[i]));
	}
	update(indis.size());
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::remove( int num,const int *id){
	if(num>m_popsize)		throw myException("the number of particles going to be deleted should less than the m_popsize@Population::remove(const int num,const int *id=0)");
		
	if(num==m_popsize){ 
		m_pop.clear();
		m_best.clear();
	}else if(id!=0){
		for(int i=0;i<num;i++){
			for(int j=0;j<m_popsize;j++){
				if(m_pop[j].get()==nullptr|| id[i]==m_pop[j]->m_id){
					m_pop.erase(m_pop.begin()+j);
					break;
				}
			}
		}
	}
	m_popsize=m_popsize-num;
	updateIndex();
	updateMaxID();
		
	findWorst();
	
	Global::msp_global->m_totalNumIndis-=num;
	m_orderList.resize(m_popsize);
}	

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::remove(const  vector<int> &id){
	int num=id.size();
	if(num>m_popsize||num<=0)		throw myException("the number of particles going to be deleted should less than the m_popsize@Population::remove(const int num,const int *id=0)");
		
	if(num==m_popsize){ 
		m_pop.clear();
		m_best.clear();
	}else{
		for(int i=0;i<num;i++){
			for(int j=0;j<m_popsize;j++){
				if(m_pop[j].get()==nullptr|| id[i]==m_pop[j]->m_id){
					m_pop.erase(m_pop.begin()+j);
					break;
				}
			}
		}
	}
	m_popsize=m_popsize-num;
	updateIndex();
	updateMaxID();
		
	findWorst();
	
	Global::msp_global->m_totalNumIndis-=num;
	m_orderList.resize(m_popsize);
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::increaseDimension(){
	Algorithm::m_numDim=GET_NUM_DIM;
	for(int i=0;i<m_popsize;i++) m_pop[i]->increaseDimension();
	for(auto &i:m_best) i->increaseDimension();
	update(0);
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::decreaseDimension(){
	Algorithm::m_numDim=GET_NUM_DIM;
	for(int i=0;i<m_popsize;i++) m_pop[i]->decreaseDimension();
	for(auto &i:m_best) i->decreaseDimension();
	update(0);
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::changeDimension(){
	if(CAST_PROBLEM_DYN->getDirDimensionChange()==true)increaseDimension();
	else decreaseDimension();
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::printPopToFile(ofstream &out){

	for(int i=0;i<m_popsize;i++){
		m_pop[i]->representative().printToFile(out);
	}
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::updateMemory(){
	for(auto &i:m_best) i->updateMemory();
	for(int j=0;j<m_popsize;j++){
		m_pop[j]->updateMemory();
	}
	findBest();
	for(unsigned i=0;i<m_bestIdx.size();i++) updateBestArchive(m_pop[m_bestIdx[i]]->representative());
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::sort(bool mode, bool representative){
	//default value for mode=true and representative=false
	vector<bool> flag;
	flag.resize(m_popsize);
	for(int i=0;i<m_popsize;i++) {
		flag[i]=0;
	}
	int index;
	for(int i=0;i<m_popsize;i++){
		int j=0;
		while(j<m_popsize&&flag[j]==1) j++;
		index=j;

		for(j=index;j<m_popsize;j++){
			if(mode&&representative){
				if(flag[j]==0&&m_pop[j]->representative()>(m_pop[index]->representative())) index=j;
			}else if(!mode && representative){
				if(flag[j]==0&&m_pop[j]->representative()<(m_pop[index]->representative())) index=j;

			}else if(mode&&!representative){
				if(flag[j]==0&&m_pop[j]->self()>(m_pop[index]->self()))	 index=j;
			}else{
				if(flag[j]==0&&m_pop[j]->self()<(m_pop[index]->self()))	 index=j;
			}

		}
		m_orderList[i]=index;
		flag[index]=1;
	}
	flag.clear();
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::rank(){
	vector<int> rank_(m_popsize,0);
	vector<int> count(m_popsize,0);
	vector<vector<int> > cset(m_popsize,vector<int>(m_popsize));
	
	for (int i = 0; i < m_popsize; i++)
	{
		m_pop[i]->m_ranking = -1;	
	}
	
	#ifdef OFEC_DEMON
		startRankThread(rank_, count, cset);
	#else
		for (int k = 0; k<m_popsize; k++)
			for (int j = 0; j<m_popsize; j++)
			{
				if (k != j)
				{
					if (*m_pop[j]>*m_pop[k]) rank_[k]++;

					if (*m_pop[k]>*m_pop[j])
					{
						cset[k][count[k]] = j;
						count[k]++;
					}
				}
			}
	#endif

	m_curRank  = 0;
	vector<int> rank2(m_popsize);
	while(1)
	{
		int stop_count = 0;
	    for(int k=0; k<m_popsize; k++)
			rank2[k] = rank_[k];

        for(int k=0; k<m_popsize; k++)
		{			
		    if(m_pop[k]->m_ranking==-1&&rank_[k]==0)
			{
				m_pop[k]->m_ranking = m_curRank;
				for(int j=0; j<count[k]; j++)
				{
				   int id =	cset[k][j];
				   rank2[id]--;
				   stop_count++;
				}
			}						
		}

	    for(int k=0; k<m_popsize; k++)
			rank_[k] = rank2[k];

		m_curRank++;
		if(stop_count==0) 
			break;
	}

}

template <typename ED, typename TypeIndi>
double Population<ED, TypeIndi>::getMaxObj(const int idx){
	double obj;
	obj=m_pop[0]->representative().obj(idx);

	for(int i=1;i<m_popsize;i++){
		if(obj<m_pop[i]->representative().obj(idx)){
			obj=m_pop[i]->representative().obj(idx);
		}
	}
	return obj;
}

template <typename ED, typename TypeIndi>
double Population<ED, TypeIndi>::getMinObj(const int idx){
	double obj;
	obj=m_pop[0]->representative().obj(idx);

	for(int i=1;i<m_popsize;i++){
		if(obj>m_pop[i]->representative().obj(idx)){
			obj=m_pop[i]->representative().obj(idx);
		}
	}
	return obj;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::checkOverCrowd(int subSize){
	if(this->m_popsize<= subSize) return;
	for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->m_flag=0;
	int num=this->m_popsize-subSize;
	int *id=new int[num];
	for(int i=0;i<num;i++){
		int index,j=0;
		while(this->m_pop[j]->m_flag!=0&& j<this->m_popsize)j++;
		index=j;

		for(int j=index+1;j<this->m_popsize;j++){
			if(this->m_pop[j]->m_flag==0&&!(this->m_pop[j]->representative()>(this->m_pop[index]->representative()))){
				index=j;
			}
		}
		id[i]=this->m_pop[index]->m_id;
		this->m_pop[index]->m_flag=1;
	}
	remove(num,id);
	delete [] id;
	id=0;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::ROOTBest(vector<int> idx){
	if(mROOT::getROOT()){		
		for(unsigned int j=0;j<idx.size();j++) mROOT::getROOT()->record(*m_pop[idx[j]]);
	}
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::ROOTBest(){
	if(mROOT::getROOT()){		
		for(auto &i:m_best) mROOT::getROOT()->record(i->self(),Global::msp_global.get());
	}
}

template <typename ED, typename TypeIndi>
const vector<unique_ptr<TypeIndi>> & Population<ED, TypeIndi>::getPop(){
	return m_pop;
}

template <typename ED, typename TypeIndi>
const int Population<ED, TypeIndi>::getPopSize()const{
	return m_popsize;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::handleReturnFlag(ReturnFlag f){
	switch (f){
	case Return_Normal:
		break;
	case Return_Change:
		updateMemory();
		break;
	case Return_Terminate:
		break;
	case Return_ChangeNextEval:
		ROOTBest();
		break;
	case Return_Change_Timelinkage:
		CAST_PROBLEM_DYN->getTriggerTimelinkage()=false;
		updateMemory();
		break;
	case Return_Loop:
		break;
	case Return_Change_Dim:
		changeDimension();
		break;
	default:
		break;
	}
}
template <typename ED, typename TypeIndi>
bool Population<ED, TypeIndi>::isConverged(){
	return false;
}
	

template <typename ED, typename TypeIndi>
TypeIndi* Population<ED, TypeIndi>::operator[](const int i){
	return m_pop[i].get();
}

template <typename ED, typename TypeIndi>
const TypeIndi* Population<ED, TypeIndi>::operator[](const int i)const{
	return m_pop[i].get();
}

template <typename ED, typename TypeIndi>
bool Population<ED, TypeIndi>::operator >(const Population<ED, TypeIndi>& p){
// comparison in terms of the number of dominated solutions of m_best
	if(m_popsize==0&&p.m_popsize>=0) return false;
	else if(m_popsize>0&&p.m_popsize==0) return true;
	
	int num1=0,num2=0;
	for(auto &i:m_best){
		for(auto& j:p.m_best){
			if(*i>*j) num1++;
			if(*i<*j) num2++;
		}
	}
	if(num1>num2) return true;
	return false;
}

template <typename ED, typename TypeIndi>
bool Population<ED, TypeIndi>::operator <(const Population<ED, TypeIndi>& p){
// comparison in terms of the number of dominated solutions of m_best
	if(m_popsize==0&&p.m_popsize>=0) return true;
	else if(m_popsize>0&&p.m_popsize==0) return false;
	int num1=0,num2=0;
	for(auto &i:m_best){
		for(auto& j:p.m_best){
			if(*i<*j) num1++;
			if(*i>*j) num2++;
		}
	}
	if(num1>num2) return true;
	return false;
}

template <typename ED, typename TypeIndi>
const vector<unique_ptr<Solution<ED>>>& Population<ED, TypeIndi>::getBest()const{
	return m_best;
}

template <typename ED, typename TypeIndi>
int Population<ED, TypeIndi>::getID(){
	return m_popID;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::setID(int id){
	m_popID=id;
}
template <typename ED, typename TypeIndi>
vector<bool> Population<ED, TypeIndi>::getFlag() const { 
	return m_flag; 
}
template <typename ED, typename TypeIndi>
int Population<ED, TypeIndi>::getEvoNum() const { 
	return m_evoNum; 
}

template <typename ED, typename TypeIndi>
int Population<ED, TypeIndi>::getCurRank(){
	return m_curRank;
}

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::resize(int size){
	if (size > m_popsize){
		for (int i = m_popsize; i < size; i++){
			m_pop.push_back(unique_ptr<TypeIndi>(new TypeIndi()));
		}
	}else if (size < m_popsize){
		m_pop.resize(size);
	}
	m_popsize = size;
}

template <typename ED, typename TypeIndi>
double Population<ED, TypeIndi>::rank(const ED & s, bool mode){
	// get ranking of a solution in terms of a ranked pop
	int rank_ = 0;
	int count_d = 0;
	int count_n = 0;
	for (int i = 0; i < m_popsize; ++i){
		if (mode){
			if (m_pop[i]->self()>s){
				if (m_pop[i]->m_ranking > rank_) rank_ = m_pop[i]->m_ranking;
				count_d++;
			}
			else if( m_pop[i]->self()<s)			count_n++;
		}else{
			if (m_pop[i]->representative()>s){
				if (m_pop[i]->m_ranking > rank_) rank_ = m_pop[i]->m_ranking;
				count_d++;
			}
			if (m_pop[i]->representative()<s) count_n++;
		}
	}
	if (count_n == m_popsize) return -1;
	if (count_n>0&&count_n + count_d == m_popsize) return rank_ + 0.5;
	if (count_n == 0 && count_d == 0) return rank_;
	return rank_ + 1;
	
}

#ifdef OFEC_DEMON
template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::startRankThread(vector<int> &rank_, vector<int> &count, vector<vector<int> > &cset){


	typedef boost::packaged_task<int> TASK;
	typedef boost::unique_future<int> FUTURE;

	int numTask = std::thread::hardware_concurrency();

	if (numTask>m_popsize) numTask = m_popsize;

	vector<TASK>  at(numTask);
	vector<FUTURE> af(numTask);

	vector<vector<int>> tsk(numTask);
	int num = (m_popsize / numTask);
	int i = 0;
	for (; i < m_popsize; ++i){
		int k = i/num;
		if (k == numTask) break;		
		tsk[k].push_back(i);
	}
	for (int j=0; i < m_popsize; ++i,++j){
		tsk[j].push_back(i);
	}
	for (int i = 0; i<numTask; i++){
		at[i] = TASK(boost::bind(Population<ED, TypeIndi>::rankThread, boost::ref(tsk[i]), boost::ref(*this), boost::ref(rank_), boost::ref(count), boost::ref(cset)));
		af[i] = at[i].get_future();
		boost::thread(boost::move(at[i]));

	}
	boost::wait_for_all(af.begin(), af.end());
	for (auto &i : af) assert(i.is_ready() && i.has_value());
}
template <typename ED, typename TypeIndi>
int Population<ED, TypeIndi>::rankThread(vector<int> &tsk, const Population<ED, TypeIndi> &p, vector<int> &rank_, vector<int> &count, vector<vector<int> >& cset){

	for (auto k:tsk){
		for (int j = 0; j<p.m_popsize; j++)
		{
			if (k != j)
			{
				if (*(p.m_pop[j])>*(p.m_pop[k])) rank_[k]++;

				if (*(p.m_pop[k])>*(p.m_pop[j]))
				{
					cset[k][count[k]] = j;
					count[k]++;
				}
			}
		}
	}
	return 0;
}
#endif

template <typename ED, typename TypeIndi>
void Population<ED, TypeIndi>::clearBestArchive(){
	m_best.clear();
}
#endif // POPULATION_H
