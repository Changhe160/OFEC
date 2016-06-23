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

#ifndef ANT_H
#define ANT_H

#include "../../Utility/definition.h"
#include "../../Problem/Combination/TSP/TravellingSalesman.h"
#include "../Individual.h"

class Ant: public Individual<CodeVInt>
{
public:
	using Individual<CodeVInt>::initialize;
	Ant();
	~Ant();
	virtual void selectNextCity_Greedy(const vector<vector<double> > &phero,double beta,double alpha=1);
	virtual void selectNextCity_Pro(const vector<vector<double> > &phero,double beta,double alpha=1); 
	//select the next city from a candidate list
	virtual void selectNextCity_Pro(const vector<vector<double> > &phero, const vector<set<int>> &candidate, double beta, double alpha = 1);

	void selectNext(Individual<CodeVInt> &parent, const vector<vector<double>> &phero, double);
	void selectNext(Individual<CodeVInt> &parent, const vector<double> &phero, double);
	void initialize(int initNode=-1);
	pair<int, int> getCurrentEdge();
	pair<int, int> getLastEdge();
	void resetData(bool flag = true);
	void initializeIteIndivl(const Solution<CodeVInt> &parent);
	void selectNextCity_GLMemory_TSP(Individual<CodeVInt> &parent, const vector<vector<double> > &phero, double);
	void selectNextCity_GLMemory_QAP(Individual<CodeVInt> &parent, const vector<vector<double> > &phero, double);
	void selectNextCity_GLMemory_MKP(Individual<CodeVInt> &parent, const vector<double> &phero, double);
protected:
	vector<int> mv_tabu; 
	int m_loc;
	vector<double> mv_accPro;
	Individual<CodeVInt> m_iteIndivl;
	void modifiedAnt(Individual<CodeVInt> &parent);
};

#endif