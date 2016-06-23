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

Thomas Stutzle. (2000). "Max¨Cmin ant system" FutureGener.Comput.Syst.
*************************************************************************/
// Created: 8 Nov 2015
// Last modified:

#ifndef MMAS_H
#define MMAS_H

#include "../Ant.h"
#include "../../Population.h"
#include "../AS/AS.h"

//for symmetric TSP

class MMAS : public AS
{
private:
	double m_pheroMax;
	double m_pheroMin;
	double m_rho;
	double m_branchFrc;
	double m_branchFactor;
	long int m_uGB;
	double m_lambda;
	int m_length; //len of candidate list

	long int m_restartFoundBest;

	Solution<CodeVInt> m_globalBest; //the globally best tour from the beginning of the trial
	Solution<CodeVInt> m_restartBest; 
	bool m_isHaveGlobalBest; // initially, there has no history best tour
	bool m_isHaveRestartBest;
	int m_impRadio;
public:
	MMAS(double alpha, double beta, double Q, int Popsize, int NC, int numDim, double coeff);
	MMAS(ParamMap &v);
	~MMAS();
	ReturnFlag run_();
	void setDefaultParameters();
	double nodeBranching();
	void checkPheromoneTrailLimits();
	void initPhero();
	void initializeSystem(double t0, bool isDefault = true);
	void updatePheromeno();
	void updatePheroMinAndMax();
	void findnearghbor();
	void pheroSmoothing(); //pheromone trail smoothing(PTS)
	double getLenOfNN();   //get the length of tour which is created based on nearest neighbor heuristic
};


#endif //MMAS_H