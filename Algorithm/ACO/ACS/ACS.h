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

Dorigo, M. (1996). "Ant Colony System: A Cooperative Learning Approach to the Traveling Salesman Problem"
IEEE Transactions on Evolutionary Computation.

*************************************************************************/
// Created: 1 Nov 2015
// Last modified:


#ifndef ANTCOLONYSYSTEM_H
#define ANTCOLONYSYSTEM_H

#include "../Ant.h"
#include "../../Population.h"

//for symmetric TSP

class ACS : public Population<CodeVInt, Ant>
{
private:
	enum StopCriterion{ MAX_GEN, MIN_COVER };
	vector<vector<double> > mvv_phero;
	double m_beta;
	double m_Q;
	int m_NC;
	double m_coeff;
	double m_t;

	int m_saveFre;
	int m_num;
	StopCriterion m_stopCriterion;
	int m_impRadio;

	Solution<CodeVInt> m_globalBest; //the globally best tour from the beginning of the trial
	bool m_isHaveGlobalBest; // initially, there has no history best tour
public:
	ACS(double beta, double Q, int Popsize, int NC, int numDim, double coeff);
	ACS(ParamMap &v);
	~ACS();
	ReturnFlag run_();
	virtual void initializeSystem();
	virtual void local_updatePheromeno(int ant_loc, bool isLastEdge = false);
	virtual void global_updatePheromeno();
	void resetAntsInfo();
	bool ifTerminating();
	double getLenOfNN();   //get the length of tour which is created based on nearest neighbor heuristic  
};


#endif //ANTCOLONYSYSTEM_H