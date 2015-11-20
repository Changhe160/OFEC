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
// Created: 21 September 2014
// Last modified: 12 Dec. 2014
#ifndef SOLUTION_H
#define SOLUTION_H

#include "global.h"
#include "encoding.h"
#include "../Utility/definition.h"
#include "../Utility/myexcept.h"
#include "../Measure/mSingleObj.h"

template<typename ED=EncodingType >
class Solution: protected ED{
public:
	Solution(int numDim,int numObj):ED(numDim,numObj){}
	Solution(const ED& rhs):ED(rhs){}
	Solution(const Solution<ED>&rhs) :ED(rhs.data()){}
	Solution():ED((Global::msp_global->mp_problem!=nullptr)?Global::msp_global->mp_problem->getNumDim():0,\
		(Global::msp_global->mp_problem!=nullptr)?Global::msp_global->mp_problem->getNumObj():0){}
	Solution & operator=(const Solution<ED>& );	
	
	ReturnFlag initialize(const bool mode=true,const int rIndex=0,const int rSize=0);
	ReturnFlag initialize( const VirtualEncoding & rhs,double radius=1.,bool mode=true);
	inline double obj(int i)const ;
	inline const vector<double> & obj(void)const;
	inline vector<double> & obj(void);
	inline ED & data();
	inline const ED & data()const;
	inline Solution<ED> & self();
	inline const Solution<ED> & self() const;
	inline int getNumDim()const ;
	double getDistance(const Solution<ED> & rhs,DistanceMode mode=DIS_EUCLIDEAN)const;
	double getDistance(const ED & rhs,DistanceMode mode=DIS_EUCLIDEAN)const;
	double getObjDistance_(const vector<double>&rhs)const{
		return this->getObjDistance(rhs);
	}
	void printToScreen();
	void printToScreen()const;
	virtual void printToFile(ofstream & out);
	void validate(SolutionValidation *mode=0);
	ReturnFlag evaluate(const bool rFlag=true);

	/// begin comparison regarding objective
	bool operator>=(const Solution<ED> &p)const{
		return this->operator>=(p.data());
	}
	bool operator<=(const Solution<ED> &p)const{
		return this->operator<=(p.data());
	}
	bool operator>(const Solution<ED> &p)const{
		return this->operator>(p.data());
	}
	bool operator<(const Solution<ED> &p)const{
		return this->operator<(p.data());
	}
	bool operator==(const Solution<ED> &p)const{
		return this->operator==(p.data());
	}
	bool operator!=(const Solution<ED> &p)const{
		return this->operator!=(p.data());
	}
	CompareResultFlag compare_(const Solution<ED> &rhs)const{
		return compare_(rhs.data());
	}
	bool operator>=(const ED &p)const; //better than or equal to p
	bool operator<=(const ED &p)const; //worse than or equal to p
	bool operator>(const ED &p)const;  //better than p
	bool operator<(const ED &p)const;  //worse than p
	bool operator==(const ED &p)const; //equal to p
	bool operator!=(const ED &p)const; //not equal to p
	inline CompareResultFlag compare_(const ED &rhs)const;

	bool isSame(const Solution &p)const;
	bool isValid();
	//// the end
	virtual void increaseDimension();
    virtual void decreaseDimension();
	virtual void updateMemory();
	Solution(Solution<ED>&& rhs);
	Solution<ED>& operator=(Solution<ED>&& rhs);
private:
	#ifdef OFEC_CONSOLE
	static boost::thread_specific_ptr<Solution<ED>> m_bestSoFar,m_worstSoFar; // one non-dominated worst/best solution
	#endif
	#ifdef OFEC_DEMON
	static unique_ptr<Solution<ED>> m_bestSoFar,m_worstSoFar;
	#endif	
public:
	static Solution<ED>& getBestSolutionSoFar(){	
		return *m_bestSoFar;		
	}
	static Solution<ED>& getWorstSolutionSoFar(){	
		return *m_worstSoFar;		
	}
	static void allocateMemoryWB(int,int);
	static void initilizeWB(const Solution<ED>&);
	static void initilizeWB(const ED&);
};


#ifdef OFEC_CONSOLE
template<typename ED >
boost::thread_specific_ptr<Solution<ED>> Solution<ED>::m_bestSoFar(nullptr);
template<typename ED >
boost::thread_specific_ptr<Solution<ED>> Solution<ED>::m_worstSoFar(nullptr); 
#endif

#ifdef OFEC_DEMON
template<typename ED >
unique_ptr<Solution<ED>> Solution<ED>::m_bestSoFar(nullptr);
template<typename ED >
unique_ptr<Solution<ED>> Solution<ED>::m_worstSoFar(nullptr);
#endif

template<typename ED>
void Solution<ED>::allocateMemoryWB(int numDim,int numObj){
	
	g_mutexStream.lock();	
	m_bestSoFar.reset(new Solution<ED>(numDim, numObj));
	m_worstSoFar.reset(new Solution<ED>(numDim, numObj));		
	g_mutexStream.unlock();
}

template<typename ED>
void Solution<ED>::initilizeWB(const Solution<ED>& s){
	*m_bestSoFar=*m_worstSoFar=s;
}
template<typename ED>
void Solution<ED>::initilizeWB(const ED& s){
	m_bestSoFar->data()=m_worstSoFar->data()=s;
}
template<typename ED>
Solution<ED> & Solution<ED>::operator=(const Solution<ED>& rhs){
	if(&rhs==this) return *this;
	if(this->m_x.size()!=rhs.m_x.size()) throw myException("m_x@Solution & operator=(const Solution & )");
	if(this->m_obj.size()!=rhs.m_obj.size()) throw myException("m_obj@Solution & operator=(const Solution &)");	
	ED::operator=(rhs);
	return *this;
}

template<typename ED>
ReturnFlag Solution<ED>::evaluate(const bool rFlag){
	ReturnFlag rf=Return_Normal;
	
	rf=Global::msp_global->mp_problem->evaluate(*this,rFlag);
	//cout<<"Evals: "<<Global::msp_global->mp_problem->getTevals()<<" "<<Global::msp_global->mp_problem->getEvaluations()<<" Flag "<<rFlag<<endl;
	//printToScreen();
	//if(Global::msp_global->mp_problem->getTevals()>=28753) getchar();

	#ifdef OFEC_DEMON 
	if (rf==Return_Error)	this->m_obj=m_worstSoFar->data().m_obj;
	#endif
	#ifdef OFEC_CONSOLE
		if(rFlag){
			if(mSingleObj::getSingleObj()!=nullptr){
				mSingleObj::getSingleObj()->record(Global::msp_global.get(),this->m_obj[0]);
			}
		}
	#endif

	if(rFlag){
		if(Global::msp_global->mp_problem->getEvaluations()-1==0)	*m_worstSoFar=*m_bestSoFar=*this;		
		else{	
			if(*this>*m_bestSoFar) *m_bestSoFar=*this;
			if(*this<*m_worstSoFar) *m_worstSoFar=*this;
		}
	}

	return rf;
}
 
template<typename ED>
ReturnFlag Solution<ED>::initialize(const bool mode,const int rIndex,const int rSize){
	Global::msp_global->mp_problem->initializeSolution(*this,rIndex,rSize);
	return evaluate(mode);
}

template<typename ED>
ReturnFlag Solution<ED>::initialize( const VirtualEncoding & rhs,double radius,bool mode){
	Global::msp_global->mp_problem->initializeSolution(rhs,*this,radius);
	return evaluate(mode);
}
template<typename ED>
double Solution<ED>::obj(int i)const {
	return this->m_obj[i];
}

template<typename ED>
const vector<double> & Solution<ED>::obj(void)const{
	return this->m_obj;
}
template<typename ED>
vector<double> & Solution<ED>::obj(void){
	return this->m_obj;
}
template<typename ED>	
int Solution<ED>::getNumDim()const {
	//users need to implement size() function for non-container-based encoding 
	return this->m_x.size();
}

template<typename ED>
double Solution<ED>::getDistance(const Solution<ED> & rhs,DistanceMode mode)const{
	return Global::msp_global->mp_problem->getDistance(*this,rhs,mode);
}
template<typename ED>
double Solution<ED>::getDistance(const ED & rhs,DistanceMode mode)const{
	return Global::msp_global->mp_problem->getDistance(*this,rhs,mode);
}

template<typename ED>
void Solution<ED>::increaseDimension(){
	unsigned oldDim=this->m_x.size();
	unsigned dim=Global::msp_global->mp_problem->getNumDim();
	this->m_x.resize(dim); //users need to implement resize() function for non-container-based encoding
	Global::msp_global->mp_problem->initializePartSolution(*this,oldDim,dim);
	evaluate(false);
}
template<typename ED>
void Solution<ED>::decreaseDimension(){
	int dim=Global::msp_global->mp_problem->getNumDim();
	this->m_x.resize(dim); //users need to implement resize() function for non-container-based encoding
	evaluate(false);
}
template<typename ED>
void Solution<ED>::printToScreen(){
	//not for general this->m_x, only for stl::container based this->m_x
	for(auto &i:this->m_x){cout<<i<<" ";}
	cout<<": ";
	for(auto &i:this->m_obj){ cout<<i<<" ";}
	cout<<endl;
}
template<typename ED>
void Solution<ED>::printToScreen()const{
	for(auto &i:this->m_x){cout<<i<<" ";}
	cout<<": ";
	for(auto &i:this->m_obj){ cout<<i<<" ";}
	cout<<endl;
}

template<typename ED>
void Solution<ED>::validate(SolutionValidation *mode){
	if(Global::msp_global->mp_problem->isValid(*this)) return;
	Global::msp_global->mp_problem->validate(*this,mode);
}

//better than or equal to p
template<typename ED>
bool Solution<ED>::operator>=(const ED &p)const{
	CompareResultFlag flag=compare_(p);
	return (flag==Compare_Better||flag==Compare_Equal);
}

//worse than or equal to p
template<typename ED>
bool Solution<ED>::operator<=(const ED &p)const{
	CompareResultFlag flag=compare_(p);
	return (flag==Compare_Worse||flag==Compare_Equal);
}


//better than p
template<typename ED>
bool Solution<ED>::operator>(const ED &p)const{
	return (Compare_Better==compare_(p));
}

template<typename ED>
bool Solution<ED>::operator<(const ED &p)const{
	return (Compare_Worse==compare_(p));
}

template<typename ED>
bool Solution<ED>::operator==(const ED &p)const{
	return (Compare_Equal==compare_(p));
}

template<typename ED>
bool Solution<ED>::operator!=(const ED &p)const{
	//note: != means non-dominated relation
	return (Compare_Non_Dominated==compare_(p));
}

template<typename ED>
bool Solution<ED>::isSame(const Solution<ED> &p)const{
	return    Global::msp_global->mp_problem->isSame(*this,p);   
}

template<typename ED>
CompareResultFlag Solution<ED>::compare_(const ED &rhs)const{
	
	if (Global::msp_global->mp_problem != nullptr){
		return Global::msp_global->mp_problem->compare(*this, rhs);
	}
	else{
		throw myException("problem instance is not created! @  Solution<ED>::compare_()");
	}	
}

template<typename ED>
ED & Solution<ED>::data(){
	return *this;
}
template<typename ED>
const ED & Solution<ED>::data() const{
	return *this;
}

template<typename ED>
Solution<ED> & Solution<ED>::self(){
	return *this;
}
template<typename ED>
const Solution<ED> & Solution<ED>::self() const{
	return *this;
}

template<typename ED>
void Solution<ED>::updateMemory(){
    evaluate(false);
}


template<typename ED>
 void Solution<ED>::printToFile(ofstream & out){
    for(auto &i:this->m_x){out<<i<<" ";}
	for(auto &i:this->m_obj){ out<<i<<" ";}
	out<<endl;
}
template<typename ED>
bool Solution<ED>::isValid(){
	return Global::msp_global->mp_problem->isValid(*this);
}

template<typename ED>
inline Solution<ED>::Solution(Solution<ED>&& rhs):ED(move(rhs)){

}

template<typename ED>
inline Solution<ED>& Solution<ED>::operator=(Solution<ED>&& rhs){
	if (this->m_x.size() != rhs.m_x.size()) throw myException("m_x@Solution & operator=( Solution && )");
	if (this->m_obj.size() != rhs.m_obj.size()) throw myException("m_obj@Solution & operator=( Solution &&)");
	ED::operator=(move(rhs));
	return *this;
}
#endif
