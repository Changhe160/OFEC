/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
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

#ifndef SLPARTICLE_H
#define SLPARTICLE_H

#include "../Particle.h"
#include "Progress.h"

class SLParticle :public Particle
{
public:
    SLParticle();
    virtual ~SLParticle();
    SLParticle(const SLParticle& other);
    SLParticle& operator=(const SLParticle& other);

    void setSelRatio();
    void nonLearnToLearn();
    void learnToNonLearn();
    void updateSelectionRatioMonitor();
    void updateSelectionRatioProg();
    int selectOperator();
    ReturnFlag move(const Solution<CodeVReal> &lbest,double w, double c1);

public:
    //bool m_flag; inherited from Particle , flag of Learning to gbest;
    unsigned int m_itersUnimpr;
    unsigned int m_updateFre;
    double m_learnRatio;

    static const int ms_numOperators=4;
	vector<Progress> mv_prog, mv_monitor;
};


#endif //SLPARTICLE_H