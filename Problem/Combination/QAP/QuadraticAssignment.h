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

#ifndef QUADRATIC_ASSIGNMENT_PROBLEM_H
#define QUADRATIC_ASSIGNMENT_PROBLEM_H

#include "../../problem.h"
#include "../../optimum.h"
#include "../../../Global/boundary.h"
#define CAST_QAP dynamic_cast<QuadraticAssignment *>(Global::msp_global->mp_problem.get())
class QuadraticAssignment :public Problem
{
protected:
	vector<vector<double> > mvv_flow;
	vector<vector<double> > mvv_distance;
	string m_fileName;
	BoundaryQAP m_searchRange;
	Optima<CodeVInt> m_globalOpt;
public:
	QuadraticAssignment(ParamMap& v);
	~QuadraticAssignment();
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true);
	bool isValid(const VirtualEncoding &s);
	void validate(VirtualEncoding &s, SolutionValidation *mode = 0) {}
	void readProblem();    //read source data from file
	bool getObjGlobalOpt(vector<double> &opt);
	bool getObjGlobalOpt(vector<vector<double>> &opt);
	const Optima<CodeVInt> & getGOpt()const;
	Optima<CodeVInt> & getGOpt();
	bool isGlobalOptKnown();
	QuadraticAssignment *getTypePtr();
	QuadraticAssignment &getTypeRef();
	bool isGlobalOptFound();
protected:
	void setObjSet();
	void initializeSolution(VirtualEncoding &result, const int idx = 0, const int maxId = 0);
	void initializeSolution(const VirtualEncoding &ref, VirtualEncoding &result, double range) {}
	void initializePartSolution(VirtualEncoding &result, int begin, int end) {}
	double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode);
	bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2);
};

inline bool QuadraticAssignment::isGlobalOptFound() {
	if (isGlobalOptKnown()) {
		if (m_globalOpt.isAllFound()) return true;
	}
	return false;
}


#endif //QUADRATIC_ASSIGNMENT_PROBLEM_H