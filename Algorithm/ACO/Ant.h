/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
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
	Ant();
	~Ant();
	virtual void selectNextCity_Greedy(const vector<vector<double> > &phero,double beta,double alpha=1);
	virtual void selectNextCity_Pro(const vector<vector<double> > &phero,double beta,double alpha=1); 
	virtual void selectNextCity_GLMemory(const Solution<CodeVInt> &parent,const vector<vector<double> > &phero,double); 
	void initialization(int initNode=-1);
	void resetData();
protected:
	vector<int> mv_tabu; 
	int m_loc;
private:
	vector<double> mv_accPro;
};

#endif