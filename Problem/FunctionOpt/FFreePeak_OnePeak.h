/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*******************************************************************************************/
// Created: 22 June 2015
// Last modified:
#ifndef FREEPEAK_ONEPEAK_H
#define FREEPEAK_ONEPEAK_H
#include "FFreePeak.h"
#include "FOnePeak.h"
class FFreePeak_OnePeak:public FreePeak{
public:
	FFreePeak_OnePeak(ParamMap &v);
	~FFreePeak_OnePeak(){};	
protected:
	void createPeaks_();
	void setHeight(const vector<double>&h);
	ReturnFlag evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	double getErr(){		
		return mc_maxHeight - Solution<CodeVReal>::getBestSolutionSoFar().data().m_obj[0];
	}
	double getGPR(){
		return m_globalOpt.getNumGOptFound()*1. / m_numGOpt;
	}
	//
	void outPeakRate();
	void simulateFun();
protected:
	int m_numGOpt;
	double mc_maxHeight = 100;
	vector<int> m_maxPeak;
	bool m_simulateFun = false;
};
#endif
