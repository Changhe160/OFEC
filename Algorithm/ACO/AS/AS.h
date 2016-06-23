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

Dorigo, M. (1996). "Ant system optimization by a colony of cooperating agents." 
IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS.

*************************************************************************/
// Created: 7 Oct 2014
// Last modified:

#ifndef ANTSYSTEM_H
#define ANTSYSTEM_H

#include "../Ant.h"
#include "../../Population.h"

//for symmetric TSP

class AS: public Population<CodeVInt,Ant>
{
protected:
	enum StopCriterion{MAX_GEN,MIN_COVER};
	vector<vector<double> > mvv_phero;
	double m_alpha,m_beta;
	double m_Q;
	int m_NC;
	double m_coeff;

	int m_saveFre;
	int m_num;
	StopCriterion m_stopCriterion;
	int m_impRadio;
public:
	AS(double alpha,double beta,double Q,int Popsize,int NC,int numDim,double coeff);
	AS(ParamMap &v);
	~AS();
	ReturnFlag run_();
	virtual void initializeSystem(double t0 = -1, bool isDefault = true);
	virtual void updatePheromeno();
	void resetAntsInfo();
	bool ifTerminating();	
};


#endif