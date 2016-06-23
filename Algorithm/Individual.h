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
// Created: 21 Aug. 2014
// Last modified:
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
template <typename, typename> class Population;
#include "../Global/solution.h"
template<typename ED>
class Individual: public Solution<ED>{
	template <typename, typename> friend class Population;
    public:
        explicit Individual():Solution<ED>(),m_id(-1),m_index(-1),m_ranking(-1),m_flag(-1), m_impr(false){}
        virtual ReturnFlag initialize(int rIdx,int rID, int rPopsize,bool mode=true);
        virtual ReturnFlag initialize(int rIdx,int rID,bool mode=true);
        virtual ReturnFlag initialize(bool mode=true);
        virtual void initialize(const Solution<ED> &x,int rIdx,int rID);
		virtual ReturnFlag initialize(const Solution<ED>  &p,double radius, int rIdx, int rID,bool mode);
		
		virtual void printToScreen();
        
        Individual(const Individual<ED> &rIndi);
		Individual(const Solution<ED> &rIndi);
        Individual<ED> & operator=(const Individual<ED> &rIndi);
       
		inline void setID(const int id);
		inline void setIdIdx(const int id, const int idx);
		inline int getId()const;
		inline void setFlag(int value);
		inline int getFlag() const;
		inline int & rank();
		/// begin comparison regarding objective
		virtual bool operator>=(const Individual<ED> &p)const; //better than or equal to p
		virtual bool operator<=(const Individual<ED> &p)const; //worse than or equal to p
		virtual bool operator>(const Individual<ED> &p)const;  //better than p
		virtual bool operator<(const Individual<ED> &p)const;  //worse than p
		virtual bool operator==(const Individual<ED> &p)const; //equal to p
		virtual bool operator!=(const Individual<ED> &p)const; //not equal to p
		virtual bool isSame(const Individual<ED> &p)const;
		virtual Solution<ED> & representative();
		virtual const Solution<ED> & representative() const;
		bool isImproved();
		bool isActive();
		void setActive(bool value);
		void setImpr(bool value);
		double & fitness();
protected:
		double m_fitness=0;               // fitness of each individual
        int m_id,m_index;
        int m_ranking=-1;
		int m_flag;				// for a certain purpose
		bool m_impr,m_active=true;			// if the indi gets improved, added 2016.4.20.
};


 template<typename ED>
ReturnFlag Individual<ED>::initialize(int rIdx,int rID, int rPopsize,bool mode){
    m_id=rID;
    m_index=rIdx; 
	return Solution<ED>::initialize(mode,rIdx,rPopsize);          
}
template<typename ED>
ReturnFlag Individual<ED>::initialize(int rIdx,int rID,bool mode){
    m_id=rID;
    m_index=rIdx;
    return Solution<ED>::initialize(mode);
}
template<typename ED>
ReturnFlag Individual<ED>::initialize(bool mode){
    return Solution<ED>::initialize(mode);
}
template<typename ED>
void Individual<ED>::initialize(const Solution<ED> &x,int rIdx,int rID){
    m_id=rID;
    m_index=rIdx;
    Solution<ED>::operator=(x);
}

template<typename ED>
ReturnFlag Individual<ED>::initialize(const Solution<ED>  &p,double radius, int rIdx, int rID,bool mode){
    m_id=rID;
    m_index=rIdx;
    return Solution<ED>::initialize(p,radius,mode);
}
template<typename ED>
Individual<ED>::Individual(const Individual<ED> &rhs):Solution<ED>(rhs), m_id(rhs.m_id),m_index(rhs.m_index),\
m_ranking(rhs.m_ranking), m_flag(rhs.m_flag), m_fitness(rhs.m_fitness), m_impr(rhs.m_impr) , m_active(rhs.m_active){

}

template<typename ED>
Individual<ED>::Individual(const Solution<ED> &rhs):Solution<ED>(rhs),m_id(-1),m_index(-1),m_ranking(-1),m_flag(-1), m_impr(0){

}

template<typename ED>
Individual<ED> & Individual<ED>::operator=(const Individual &rhs){
    if(this== &rhs) return *this;
	Solution<ED>::operator=(rhs);
    m_id=rhs.m_id;
    m_index=rhs.m_index;
    m_ranking=rhs.m_ranking;
	m_flag=rhs.m_flag;
	m_fitness = rhs.m_fitness;
	m_impr = rhs.m_impr;
	m_active = rhs.m_active;
    return *this;
}

 template<typename ED>
 void Individual<ED>::printToScreen(){
	cout<<"ID: "<<m_id<<endl;
	cout<<"self:  ";
	Solution<ED>printToScreen();
}

template<typename ED>
void Individual<ED>::setIdIdx(const int id, const int idx){
	m_id=id; m_index=idx;
}

template<typename ED>
void Individual<ED>::setID(const int id){
	m_id=id;
}
template<typename ED>
int Individual<ED>::getId()const {	
	return m_id; 
}
template<typename ED>
void Individual<ED>::setFlag(int value) { 
	m_flag=value; 
}
template<typename ED>
int Individual<ED>::getFlag() const { 
	return m_flag; 
}

template<typename ED>
int & Individual<ED>::rank(){
	return m_ranking;
}

//better than or equal to p
template<typename ED>
bool Individual<ED>::operator>=(const Individual<ED> &p)const{
	CompareResultFlag flag=representative().compare_(p.representative());
	return (flag==Compare_Better||flag==Compare_Equal);
}

//worse than or equal to p
template<typename ED>
bool Individual<ED>::operator<=(const Individual<ED> &p)const{
	CompareResultFlag flag=representative().compare_(p.representative());
	return (flag==Compare_Worse||flag==Compare_Equal);
}


//better than p
template<typename ED>
bool Individual<ED>::operator>(const Individual<ED> &p)const{
	return (representative().compare_(p.representative())==Compare_Better);
}

template<typename ED>
bool Individual<ED>::operator<(const Individual<ED> &p)const{
	return (Compare_Worse==representative().compare_(p.representative()));
}

template<typename ED>
bool Individual<ED>::operator==(const Individual<ED> &p)const{
	return (Compare_Equal==representative().compare_(p.representative()));
}

template<typename ED>
bool Individual<ED>::operator!=(const Individual<ED> &p)const{
	//note: != means non-dominated relation
	return (Compare_Non_Dominated==representative().compare_(p.representative()));
}

template<typename ED>
bool Individual<ED>::isSame(const Individual<ED> &p)const{
	/*bool flag=m_self.isSame(p.m_self);
	if(flag==false) return flag;
	if(&representative()!=&m_self){
		flag=representative().isSame(p.representative());
	}
	return flag;*/
	return representative().isSame(p.representative());
}

template<typename ED>
 Solution<ED> & Individual<ED>::representative(){
    return *this;
}

template<typename ED>
const Solution<ED> & Individual<ED>::representative() const{
	return *this;
};
template<typename ED>
bool Individual<ED>::isImproved() {
	return m_impr;
}
template<typename ED>
bool Individual<ED>::isActive() {
	return m_active;
}

template<typename ED>
void Individual<ED>::setActive(bool val) {
	m_active=val;
}

template<typename ED>
double & Individual<ED>::fitness() {
	return m_fitness;
}

template<typename ED>
void Individual<ED>::setImpr(bool value) {
	m_impr = value;
}
#endif // INDIVIDUAL_H
