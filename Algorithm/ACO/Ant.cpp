/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Oct 2014
// Last modified:

#include "Ant.h"
#include "../../Global/global.h"

Ant::Ant() :Individual(), m_iteIndivl()
{
	mv_tabu.resize(getNumDim(),0);
	m_loc=0;
	m_flag=false;
	mv_accPro.resize(getNumDim());
}

Ant::~Ant()
{
	mv_tabu.clear();
}

void Ant::selectNextCity_Greedy(const vector<vector<double> > &phero,double beta,double alpha)
{
	double max=-0xfffffff;
	double temp;
	int i,pos=-1;
	for(i=0;i<getNumDim();i++)
	{
		if(mv_tabu[i]==0)
		{
			temp=pow(phero[(m_x[m_loc])][i],alpha)*pow((1./dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->getCost()[(m_x[m_loc])][i]),beta);   
			if(max<temp) 
			{
				max=temp;
				pos=i;
			}
		}
	}
	if(pos==-1)
		throw myException("selectNextCity_Greedy() error in Ant.cpp");
	++m_loc;
	m_x[m_loc]=pos;
	mv_tabu[pos]=1;
}

void Ant::selectNextCity_Pro(const vector<vector<double> > &phero,double beta,double alpha){
	
	int i,curNode=m_x[m_loc],numDim=getNumDim(),firstZero=-1;

	for(i=0;i<numDim;i++){
		if(mv_tabu[i]==0){
			mv_accPro[i]=pow(phero[curNode][i],alpha)*pow(1./dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->getCost()[curNode][i],beta);
			if(firstZero==-1) firstZero=i;
		}else{
			mv_accPro[i]=0;
		}
		if(i>0) mv_accPro[i]+=mv_accPro[i-1];
	}
	++m_loc;
	if(mv_accPro[numDim-1]>0){
		double p=Global::msp_global->mp_uniformAlg->Next()*mv_accPro[numDim-1];
		//vector<double>::iterator it= find_if(mv_accPro.begin(),mv_accPro.end(),[&](const double &i){return p<=i;});
		vector<double>::iterator it= lower_bound(mv_accPro.begin(),mv_accPro.end(),p);
		m_x[m_loc]=int (it-mv_accPro.begin());
	}else{
		m_x[m_loc]=firstZero;
	}
	mv_tabu[m_x[m_loc]]=1;
}

void Ant::selectNextCity_Pro(const vector<vector<double> > &phero, const vector<set<int>> &candidate, double beta, double alpha)
{
	int i, curNode = m_x[m_loc], numDim = getNumDim(), firstZero = -1;

	for (i = 0; i<numDim; i++){
		if (mv_tabu[i] == 0 && candidate[curNode].find(i) != candidate[curNode].end()){
			mv_accPro[i] = pow(phero[curNode][i], alpha)*pow(1. / dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->getCost()[curNode][i], beta);
			if (firstZero == -1) firstZero = i;
		}
		else{
			mv_accPro[i] = 0;
		}
		if (i>0) mv_accPro[i] += mv_accPro[i - 1];
	}
	if (mv_accPro[numDim - 1]>0){
		++m_loc;
		double p = Global::msp_global->mp_uniformAlg->Next()*mv_accPro[numDim - 1];
		//vector<double>::iterator it= find_if(mv_accPro.begin(),mv_accPro.end(),[&](const double &i){return p<=i;});
		vector<double>::iterator it = lower_bound(mv_accPro.begin(), mv_accPro.end(), p);
		m_x[m_loc] = int(it - mv_accPro.begin());
		mv_tabu[m_x[m_loc]] = 1;
	}
	else{
		if (firstZero != -1)
		{
			++m_loc;
			m_x[m_loc] = firstZero;
			mv_tabu[m_x[m_loc]] = 1;
		}
		else
			selectNextCity_Greedy(phero,beta,alpha);
	}
}

void Ant::initialize(int initNode)
{
	if(initNode!=-1) m_x[m_loc]=initNode;
	else m_x[m_loc]=int((getNumDim()-1)*Global::msp_global->mp_uniformAlg->Next());
	mv_tabu[int(m_x[m_loc])]=1;
}

pair<int, int> Ant::getCurrentEdge()
{
	if (m_loc >= 1)
		return make_pair(int(m_x[m_loc - 1]), int(m_x[m_loc]));
	else
		return make_pair(-1, -1);
}

pair<int, int> Ant::getLastEdge()
{
	return make_pair(int(m_x[m_loc]), int(m_x[0]));
}

void Ant::resetData(bool flag)
{
	for(int i=0;i<getNumDim();i++)
		mv_tabu[i]=0;
	m_loc=0;
	if (flag)
		mv_tabu[(m_x[m_loc])]=1;
	m_flag=false;  // to be evaluated
}

void Ant::initializeIteIndivl(const Solution<CodeVInt> &parent)
{
	for (int i = 0; i < getNumDim(); i++)
	{
		m_iteIndivl.data().m_x[i] = parent.data().m_x[i];
	}
	m_iteIndivl.data().m_obj = parent.data().m_obj;
}

void Ant::modifiedAnt(Individual<CodeVInt> &parent)
{
	if (m_iteIndivl.data().m_x[m_loc] != m_x[m_loc])
	{
		double obj = m_iteIndivl.data().m_obj[0];
		TravellingSalesman *ptr = dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get());
		int pos = -1;
		for (int i = m_loc + 1; i < getNumDim(); i++)
		{
			if (m_iteIndivl.data().m_x[i] == m_x[m_loc])
			{
				pos = i;
				break;
			}
		}
		obj = obj - ptr->getCost()[m_iteIndivl.data().m_x[pos]][m_iteIndivl.data().m_x[pos - 1]] - ptr->getCost()[m_iteIndivl.data().m_x[pos]][m_iteIndivl.data().m_x[(pos + 1) % getNumDim()]]\
			- ptr->getCost()[m_iteIndivl.data().m_x[m_loc]][m_iteIndivl.data().m_x[m_loc - 1]];
		obj = obj + ptr->getCost()[m_iteIndivl.data().m_x[pos]][m_iteIndivl.data().m_x[m_loc - 1]] + ptr->getCost()[m_iteIndivl.data().m_x[pos]][m_iteIndivl.data().m_x[m_loc]]\
			+ ptr->getCost()[m_iteIndivl.data().m_x[pos - 1]][m_iteIndivl.data().m_x[(pos + 1) % getNumDim()]];
		m_iteIndivl.data().m_obj[0] = obj;
		for (int i = pos - 1; i >= m_loc; i--)
		{
			m_iteIndivl.data().m_x[i + 1] = m_iteIndivl.data().m_x[i];
		}
		m_iteIndivl.data().m_x[m_loc] = m_x[m_loc];
		if ((ptr->getOptType() == MIN_OPT && obj < parent.data().m_obj[0]) || (ptr->getOptType() == MAX_OPT && obj > parent.data().m_obj[0]))
		{
			parent.data().m_obj[0] = obj;
			for (int i = 0; i < getNumDim(); i++)
				parent.data().m_x[i] = m_iteIndivl.data().m_x[i];
			parent.setFlag(true);
		}
	}
}

void Ant::selectNextCity_GLMemory_TSP(Individual<CodeVInt> &parent, const vector<vector<double> > &phero, double xp){
	
	double p=Global::msp_global->mp_uniformAlg->Next();
	if(p<xp){
		pair<int,int> pa=dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->getNextCity(parent,m_x[m_loc]);
		int node=mv_tabu[pa.first]<mv_tabu[pa.second]?pa.first:pa.second;
		
		if(mv_tabu[node]==0){
			++m_loc;
			m_x[m_loc]=node;
			mv_tabu[node]=1;
			if(!m_flag&&parent.data().m_x[m_loc]!=node) m_flag=true;
			modifiedAnt(parent);
			return;
		}
	}	

	int i,curNode=m_x[m_loc],numDim=getNumDim(),firstZero=-1;
	/*for(i=0;i<getNumDim();i++){
		if(i!=curNode&&dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->getCost()[curNode][i]==0){
			++m_loc;
			m_x[m_loc]=i;
			if(m_x[m_loc]!=parent.m_x[m_loc]) m_flag=true;
			return;
		}
	}
	for(i=0;i<getNumDim();i++){
		mv_accPro[i]=0;
		if(mv_tabu[i]==0){
			sum+=phero[curNode][i];
		}
	}
	if(sum>0){
		for(i=0;i<getNumDim();i++){
			if(mv_tabu[i]==0){
				mv_accPro[i]=phero[curNode][i];
				mv_accPro[i]=mv_accPro[i]/sum;
			}
			if(i>0) mv_accPro[i]+=mv_accPro[i-1];
		}
		double p=Global::msp_global->mp_uniformAlg->Next();
		int pos;
		for(i=0;i<getNumDim();i++){
			if(p<=mv_accPro[i]){
				pos=i;
				break;
			}
		}
		++m_loc;	
		m_x[m_loc]=pos;

	}else{
		for(i=0;i<getNumDim();i++){
			if(mv_tabu[i]==0){
				++m_loc;
				m_x[m_loc]=i;
				break;
			}
		}
	}
	if(m_x[m_loc]!=parent.m_x[m_loc]) m_flag=true;
	mv_tabu[m_x[m_loc]]=1;*/

	
	for(i=0;i<numDim;i++){
		if(mv_tabu[i]==0){
			mv_accPro[i]=phero[curNode][i];
			if(firstZero==-1) firstZero=i;
		}else{
			mv_accPro[i]=0;
		}
		if(i>0) mv_accPro[i]+=mv_accPro[i-1];
	}
	++m_loc;
	if(mv_accPro[numDim-1]>0){
		double p=Global::msp_global->mp_uniformAlg->Next()*mv_accPro[numDim-1];
		//vector<double>::iterator it= find_if(mv_accPro.begin(),mv_accPro.end(),[&](const double &i){return p<=i;});
		vector<double>::iterator it= lower_bound(mv_accPro.begin(),mv_accPro.end(),p);
		m_x[m_loc]=int (it-mv_accPro.begin());
	}else{
		m_x[m_loc]=firstZero;
	}
	if(!m_flag&&m_x[m_loc]!=parent.data().m_x[m_loc]) m_flag=true;
	modifiedAnt(parent);
	mv_tabu[m_x[m_loc]]=1;
}


void Ant::selectNextCity_GLMemory_QAP(Individual<CodeVInt> &parent, const vector<vector<double> > &phero, double xp) {
	double p = Global::msp_global->mp_uniformAlg->Next();
	if (p<xp) {
		int node = parent.data().m_x[m_loc + 1];

		if (mv_tabu[node] == 0) {
			++m_loc;
			m_x[m_loc] = node;
			mv_tabu[node] = 1;
			if (!m_flag&&parent.data().m_x[m_loc] != node) m_flag = true;
			return;
		}
	}

	int i, numDim = getNumDim(), firstZero = -1;

	for (i = 0; i<numDim; i++) {
		if (mv_tabu[i] == 0) {
			mv_accPro[i] = phero[m_loc + 1][i];
			if (firstZero == -1) firstZero = i;
		}
		else {
			mv_accPro[i] = 0;
		}
		if (i>0) mv_accPro[i] += mv_accPro[i - 1];
	}
	++m_loc;
	if (mv_accPro[numDim - 1]>0) {
		double p = Global::msp_global->mp_uniformAlg->Next()*mv_accPro[numDim - 1];
		vector<double>::iterator it = lower_bound(mv_accPro.begin(), mv_accPro.end(), p);
		m_x[m_loc] = int(it - mv_accPro.begin());
	}
	else {
		m_x[m_loc] = firstZero;
	}
	if (!m_flag&&m_x[m_loc] != parent.data().m_x[m_loc]) m_flag = true;
	mv_tabu[m_x[m_loc]] = 1;
}

void Ant::selectNextCity_GLMemory_MKP(Individual<CodeVInt> &parent, const vector<double>&phero, double xp) {
	double p = Global::msp_global->mp_uniformAlg->Next();
	if (p<xp) {
		int node = parent.data().m_x[m_loc + 1];
		++m_loc;
		m_x[m_loc] = node;
		if (!m_flag&&parent.data().m_x[m_loc] != node) m_flag = true;
		return;
	}
	else {
		++m_loc;
		double p = Global::msp_global->mp_uniformAlg->Next();
		if (p<phero[m_loc])
			m_x[m_loc] = 1;
		else 
			m_x[m_loc] = 0;
		if (!m_flag&&m_x[m_loc] != parent.data().m_x[m_loc]) m_flag = true;
	}
}

void Ant::selectNext(Individual<CodeVInt> &parent, const vector<vector<double>> &phero, double xp) {
	if (Global::msp_global->mp_problem->isProTag(TSP)) {
		selectNextCity_GLMemory_TSP(parent, phero, xp);
	}
	else if (Global::msp_global->mp_problem->isProTag(QAP)) {
		selectNextCity_GLMemory_QAP(parent, phero, xp);
	}	
}

void Ant::selectNext(Individual<CodeVInt> &parent, const vector<double> &phero, double xp) {
	if (Global::msp_global->mp_problem->isProTag(MKP)) {
		selectNextCity_GLMemory_MKP(parent, phero, xp);
	}
}

