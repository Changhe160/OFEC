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
*************************************************************************/
// Created: 18 Apr 2016
// Last modified:
//注： 该问题被转换为最小化问题, 保存的最好结果没有精确检测是否是合法解

#ifndef MULTIDIMENSIONAL_KNAPSACK_PROBLEM_H
#define MULTIDIMENSIONAL_KNAPSACK_PROBLEM_H

#include "../../problem.h"
#include "../../optimum.h"
#include "../../../Global/boundary.h"
#define CAST_MKP dynamic_cast<MultidimensionalKnapsack *>(Global::msp_global->mp_problem.get())
class MultidimensionalKnapsack :public Problem
{
protected:
	vector<double> mv_p;
	vector<vector<double> > mvv_r;
	vector<double> mv_b;
	string m_fileName;
	Optima<CodeVInt> m_globalOpt;
	int m_m;
	double m_maxP;
public:
	MultidimensionalKnapsack(ParamMap& v);
	~MultidimensionalKnapsack();
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	bool isValid(const VirtualEncoding &s);
	void validate(VirtualEncoding &s, SolutionValidation *mode = 0) {}
	void readProblem();    //read source data from file
	int invalidConstrainNum(VirtualEncoding &s_);
	bool getObjGlobalOpt(vector<double> &opt);
	bool getObjGlobalOpt(vector<vector<double>> &opt);
	const Optima<CodeVInt> & getGOpt()const;
	Optima<CodeVInt> & getGOpt();
	bool isGlobalOptKnown();
	MultidimensionalKnapsack *getTypePtr();
	MultidimensionalKnapsack &getTypeRef();
	bool isGlobalOptFound();
protected:
	void setObjSet();
	void initializeSolution(VirtualEncoding &result, const int idx = 0, const int maxId = 0);
	void initializeSolution(const VirtualEncoding &ref, VirtualEncoding &result, double range) {}
	void initializePartSolution(VirtualEncoding &result, int begin, int end) {}
	double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode);
	bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2);
};

inline bool MultidimensionalKnapsack::isGlobalOptFound() {
	if (isGlobalOptKnown()) {
		if (m_globalOpt.isAllFound()) return true;
	}
	return false;
}

#endif //MULTIDIMENSIONAL_KNAPSACK_PROBLEM_H
