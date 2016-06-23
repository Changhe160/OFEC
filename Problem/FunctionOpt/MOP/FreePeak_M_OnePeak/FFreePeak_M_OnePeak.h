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
// Created: 11 July 2015
// Last modified:
#ifndef FREEPEAK_M_ONEPEAK2_H
#define FREEPEAK_M_ONEPEAK2_H
#include "../../FFreePeak.h"

class FFreePeak_M_OnePeak :public FreePeak{
	enum CASE{ Jump = 1, Web,Countable};
	struct POS{
		CodeVReal sol;
		int order=0,pfr=0;
		bool isboundary=false,ispeak=false;
		POS(const CodeVReal &s) :sol(s){};
		POS(CodeVReal &&s) :sol(move(s)){};
	};
public:
	FFreePeak_M_OnePeak(ParamMap &v);
	ReturnFlag evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	bool isValid(const VirtualEncoding  &s);
	const vector<pair<double, double>> & PF(){ return m_pfr; }
	const vector<pair<double, double>> & PF()const { return m_pfr; }
	const vector<int>& getParetoRegion();
	const vector<bool> &getisPSR(){ return m_isPR; }
	int getType() { return m_case; }
protected:
	void createPeaks_();
	void computePFR();
	
	void configAttraction();
	void computePS();
	void getSegmentPS(int oidx, int sh1, const POS &o1, const POS &o2, vector<pair<POS, POS>> &mseg);
	static void getPSBetween(const POS &x1, POS &x2, list<POS>& ps, FFreePeak_M_OnePeak&, int maxNum,int oidx, int layer, int region, bool isPeak, bool included);
	static int PSThread(int oidx, int shape, POS &o, list<POS>& ps, vector<list<POS>::iterator> &tsk, FFreePeak_M_OnePeak&);

protected:
	int m_numParetoRegion=1;
	CASE m_case;
	vector<int> m_objShape;				// shape of each objective, the peaks of the same objective have the same shape in the Web and Jump cases
	vector<pair<double, double>>m_pfr;	//ranges of objectives on pareto front
	vector<vector<int>> m_ipsr;	//idexes of pareto solution regions, each element indicates PS regions which map to the PF
	vector<int> m_ipsr1;
	bool m_attraction=false;
	vector<bool> m_isPR;	// flags of pareto regions
	
	list<POS> m_ps;	//pareto set
	vector<vector<vector<POS>>> m_psr;	//pareto solution regions
};

#endif