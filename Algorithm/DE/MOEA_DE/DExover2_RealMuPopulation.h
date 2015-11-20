/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*  See the details of MOEA/D-DE in the following paper
*  H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization 
*  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
*  University of Essex, 2007
*************************************************************************/
// Created: 30 December 2014
// Last modified:

#ifndef DEXOVER2_REALMUPOPULATION_H
#define DEXOVER2_REALMUPOPULATION_H

#include "../DEPopulation.h"
#include "../DEIndividual.h"

class DExover2_RealMuPopulation : public DEPopulation<CodeVReal,DEIndividual>
{
public:
	DExover2_RealMuPopulation(int popsize,bool mode=true);
	DExover2_RealMuPopulation();
	~DExover2_RealMuPopulation(){}
	void cross_mutate(const vector<int> &index, DEIndividual &child);     

protected:
	double m_r;
	double m_etam;
};


#endif //DEXOVER2_REALMUPOPULATION_H