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
*************************************************************************/
// Created: 17 July 2015
// Last modified:
#ifndef FREEPEAK_D_M_ONEPEAK_H
#define FREEPEAK_D_M_ONEPEAK_H
#include "../FunctionOpt/FFreePeak.h"
class FreePeak_D_M_OnePeak :public FreePeak{
public:
	enum TYPE{ CT_Web = 1, CT_Rotation, CT_Jump, CT_Basin, CT_NumObj, CT_NumDim };
	FreePeak_D_M_OnePeak(ParamMap &v);
	ReturnFlag evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	void copyChanges(const Problem * op, const vector<int> *cd = nullptr, const vector<int> *co = nullptr);
	const vector<int>& getParetoRegion(){ return m_ipsr; }
	const vector<bool> &getisPSR(){ return m_isPR; }
	int getChangeInteval() { return m_changeInterval; }
protected:
	void createPeaks_();	
	void parameterSetting(Problem * rP);
	void jump();
	void rotate();
	void web();
	void varyBasin();
	void change();
	void varyNumDim();	//remove random dimensions when dim decreases,assume # of boxes does not change
	void varyNumObj();	//
	void varyNoBoxes(); //
protected:
	vector<double> m_radius;	// all peaks are located at the surface of a hyberball with radius of m_radius
	vector<int> m_ipsr; //the idex of the region where POF sits
	TYPE m_type;	// change type
	double m_shiftSeverity=0.1,m_heightSeverity=0.1;
	int m_centerPeak, m_onCirclePeak;
	int m_changeInterval;
	struct BI{
		size_t cutdim;		// the dimension to be changed;
		double low, high; //range to be changed
		size_t idx;		//index of the division point to be changed
	};	//information of basin to be changed
	vector<BI> m_basin;
	int m_counter=0;    // the current number of changes
	int m_initNumDim,m_initNumObj;
	double mc_maxHeight = 100;
	void update_psr();
	vector<bool> m_isPR;	// flags of pareto regions
	bool m_tl=false;		//time-linkage feature
	vector<bool> m_fcb;		//regions to change
};
#endif