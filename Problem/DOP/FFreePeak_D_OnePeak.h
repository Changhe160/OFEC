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
#ifndef FREEPEAK_D_ONEPEAK_H
#define FREEPEAK_D_ONEPEAK_H
#include "../FunctionOpt/FFreePeak.h"
#define CAST_PROBLEM_DYN_ONEPEAK dynamic_cast<FFreePeak_D_OnePeak*>(Global::msp_global->mp_problem.get())

class FFreePeak_D_OnePeak final :public FreePeak{
protected:
	void changeHeight(int &numChangeGOpt);
	void changeLocation();
	void changeNumBoxes();
	void addNoise(vector<double>& x, int boxIdx, int tidx, int gidx);
	void change();
	void createPeaks_();
	void setupGOpt();
	void updateChangingPeaksIdx(int &numChangeGOpt);
	void updateChangingBasinIdx(int &numChangeGOpt);
	void changeBasin(const int i);
	int nearestPeak(const CodeVReal &s);
	void updateMaxHeight();
	void changeShape();
public:
	FFreePeak_D_OnePeak(ParamMap &v);
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	bool predictChange(const int evalsMore);
	double getErr();
	double getRobustErr();
	double getGPR();
	double getPR();
	int getPT();
	bool isAllGOptTraced();
	void setShapeSet(vector<int> *sset=nullptr);
	int getChangeInteval() { return m_changeInterval; }
	int getChangeCounter() {return m_counter;	}
protected:
	 
	bool m_noise = false;			//whether environments are noisy
	double m_detect;			//[0,1];
	double m_lambada = 0.;
	enum Feature1{FT_timelink=1,FT_detect,FT_default};
	Feature1 m_feature1= FT_default;					//one of the above three features
	union{
		bool m_basinChange;			// whether basion boundaries change
		bool m_predic;			//whether changes is predictable
	}uf1; // mutual exclusive features
	enum Feature2{ FT_basin=1, FT_predic, FT_default2};
	Feature2 m_feature2= FT_default2;					//one of the above two features

	double m_shiftSeverity=0.01;	//[0-1];
	const double mc_maxHeightSeverity = 7;
	int m_changeInterval=5000;		// changes occur every m_changeInterval objective evaluations
	const double mc_maxHeight = 100,mc_minHeight=0;
	int m_step = 2;					// for changes of the number of peaks
	int m_numGOpt=1;
	vector<int> m_changingPeak;
	vector<bool> m_isFound;
	int m_counter=0;
	bool m_numPeakChange=false;
	struct BI{
		size_t cutdim;		// the dimension to be changed;
		double low, high; //range to be changed
		size_t idx;		//index of the division point to be changed
	};	//information of basin to be changed
	vector<BI> m_basin;
	int m_initNumPeak;
	vector<bool> m_tl;	// time-linkage flag
	int m_timewindow=1;	//for ROOT
	int m_robustPeak;

	//for evalutation
	map<int,double> m_nearSolut;

	double m_maxHeight;
	double m_noiseSeverity = 0.01;

	//for shape change
	vector<int> m_shapeSet;
};

#endif