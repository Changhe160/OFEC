#ifndef SBXREALMUPOPULATION_H
#define SBXREALMUPOPULATION_H

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
*************************************************************************/
// Created: 11 Jan 2015
// Last modified:

#include "../GAPopulation.h"
#include "../GAIndividual.h"

class SBX_RealMuPopulation : public GAPopulation<CodeVReal,GAIndividual<CodeVReal>>
{
public:
	SBX_RealMuPopulation(int popsize,bool mode=true);
	SBX_RealMuPopulation();
	~SBX_RealMuPopulation(){}
	void cross_mutate(const vector<int> &index, GAIndividual<CodeVReal> *child1, GAIndividual<CodeVReal> *child2);
	void setCrossXP(double cr) { m_cr=cr; }
	void setMutationP(double mr) { m_mr=mr; }
	void setEta(double ceta, double meta) { m_ceta=ceta; m_meta=meta; }
private:
	double m_cr;
	double m_ceta;
	double m_mr;
	double m_meta;
};

void SimulatedBinaryCrossover(GAIndividual<CodeVReal> *child1, GAIndividual<CodeVReal> *child2, const GAIndividual<CodeVReal> &parent1, const GAIndividual<CodeVReal> &parent2,  double cr, double ceta);
void PolynomialMutation(GAIndividual<CodeVReal> *indv, double mr, double meta);
#endif //SBXREALMUPOPULATION_H