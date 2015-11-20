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
// Last modified: 19/05/2013
#ifndef GROUP_H
#define GROUP_H

#include "../TypeList/TypeManip.h"
#include "../../Algorithm/Individual.h"
template<typename,typename> class Cluster;

/*TypeIndi must be a class inheriated from Indiviudal*/
template <typename ED,typename TypeIndi=Individual<Solution<ED>>>
class Group{
	template<typename,typename> friend class Cluster;
private:
	Solution<ED> m_center;
	vector<Solution<ED>> m_best;
	vector<TypeIndi *> m_member;
	int m_ID;
	double m_radius;
public:
	//Group(int runId,const int num);
	Group(const Group &g );
	Group();
	~Group(){
		m_member.clear();
	};
	void initialize( unique_ptr<TypeIndi>&  p,const int id);
	void initialize(const int num,const int id);
	void initialize( Group &g );
	bool operator ==(const Group &g);
	void merge(const Group &g);
	Group & operator=(const Group &g);
	void calculateRadius();
	int getSize()const{
		return m_member.size();
	}
	const vector<Solution<ED>>& getBest()const{
		return m_best;
	}
	TypeIndi* operator[](const int i){
		return m_member[i];
	}
	const TypeIndi* operator[](const int i)const{
		return m_member[i];
	}

	const Solution<ED> & getCenter()const{
		return m_center;
	}
	void updateBestArchive(const Solution<ED> & chr);
};
// implementation of Group class
template <typename ED,typename TypeIndi>
Group<ED,TypeIndi>::Group():m_center(),m_best(0),m_member(0),m_ID(-1),m_radius(-1){
	bool check[Loki::SuperSubclass<Individual<ED>,TypeIndi>::value]; // check if TypeIndi is derived from Individual or is type Solution
}

template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::initialize( unique_ptr<TypeIndi> &  p,const int id){
	m_member.resize(1);
	m_member[0]=p.release(); 
	m_center=m_member[0]->representative();
	m_ID=id;
	m_radius=0;
	m_best.clear();
	m_best.push_back(m_member[0]->representative());
}
template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::initialize(const int num,const int id){
	m_member.resize(num);
	m_ID=id;
}
template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::initialize( Group &g ){

	m_member.resize(g.m_member.size());
	for(unsigned i=0;i<m_member.size();i++){
		m_member[i]=g.m_member[i];
	}
	m_center=g.m_center;
	m_ID=g.m_ID;
	m_radius=g.m_radius;
	m_best=g.m_best;
}
template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::merge( const Group &g){

	for(unsigned i=0;i<g.m_member.size();i++){
		m_member.push_back( g.m_member[i]);
	}
	for(auto& i:g.m_best){
		updateBestArchive(i);
	}
	calculateRadius();
}
template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::updateBestArchive(const Solution<ED> & chr){
		bool first=true;
		// check dominated case
		for(auto i=m_best.begin();i!=m_best.end();i++){
			if(first && *i<chr){
				*i=chr;
				first=false;
			}else if(!first && *i<chr){
				m_best.erase(i);
				i--;
			}
		}
		if(!first) return;
		//check equal case	
		for(auto i=m_best.begin();i!=m_best.end();i++){
			if(*i==chr&&!(i->isSame(chr))){
				m_best.push_back(chr);
				return;
			}
		}
		//check non-dominated case	
		for(auto i=m_best.begin();i!=m_best.end();i++){
			if(!(*i!=chr)) return;
		}
		m_best.push_back(chr);
		
	}
template <typename ED,typename TypeIndi>
void Group<ED,TypeIndi>::calculateRadius(){

	m_radius=0;
	if(m_member.size()<2) return;

	if(Global::msp_global->mp_problem->isProTag(CONT)){
					
		for( int i=0;i<GET_NUM_DIM;i++){
			double x=0.;
			for(unsigned j=0;j<m_member.size();j++){
				Solution<ED> & chr=m_member[j]->representative();
				x=chr.data()[i]+x;
			}
			m_center.data()[i]=x/m_member.size();
		}
	}

	for(unsigned j=0;j<m_member.size();j++)
		m_radius+=m_member[j]->representative().getDistance(m_center);
	m_radius=m_radius/m_member.size();

	m_center.evaluate(false);
	updateBestArchive(m_center);
}
template <typename ED,typename TypeIndi>
bool Group<ED,TypeIndi>::operator ==(const Group &g){

	if(m_member.size()!=g.m_member.size()) return false;
	for(int i=0;i<m_member.size();i++){
		int j=0;
		while(j<m_member.size()&&m_member[i].m_id!=g.m_member[j].m_id)j++;
		if(j==m_member.size()) return false;
	}
	return true;
}
template <typename ED,typename TypeIndi>
Group<ED,TypeIndi>& Group<ED,TypeIndi>::operator =( const Group &g){

	if(this==&g) return *this;
	m_member=g.m_member;
	m_center=g.m_center;
	m_radius=g.m_radius;
	m_ID=g.m_ID;
	m_best=g.m_best;
	return *this;
}
template <typename ED,typename TypeIndi>
Group<ED,TypeIndi>::Group( const Group &g ):m_center(g.m_center),m_best(g.m_best),m_member(g.m_member){	
	m_ID=g.m_ID;
	m_radius=g.m_radius;
}
#endif // GROUP_H
