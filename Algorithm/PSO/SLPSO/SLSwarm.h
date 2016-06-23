/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 19 Jan 2015
// Last modified:

/*C. Li,S. Yang, and T. T. Nguyen. A self-learning particle swarm optimizer for 
global optimization problems.IEEE Transactions on Systems, Man, and Cybernetics 
Part B: Cybernetics*/

#ifndef SLPSO_H
#define SLPSO_H

#include "SLParticle.h"
#include "../Swarm.h"

class SLPSO :public Swarm<CodeVReal,SLParticle>
{

public:
    int m_numLearnToGbest;
    static int ms_updateFre;
    static double ms_learnRatio;
    static float ms_ratioLearnToGbest;
public:
	SLPSO(ParamMap &v);
	~SLPSO(){}
    void updateLearnToGbest();
	ReturnFlag evolve();
    ReturnFlag run_();
    void setParameters();
    void updateParameters();
    int getNumLearnTogbest();
    void setNumLearning();
    void calculateNumLearning(const int sfes);
    void MRandToBest( int num);
};

#endif //SLPSO_H