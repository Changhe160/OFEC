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


#include "MultidimensionalKnapsack.h"
#include "../../../Global/global.h"
#include<string>

using namespace std;

MultidimensionalKnapsack::MultidimensionalKnapsack(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1), m_globalOpt(m_numDim, m_numObj)
{
	v[param_numObj] = 1;
	m_fileName = (string)v[param_dataFile1];
	readProblem();

	int mode = int(v[param_populationInitialMethod]);
	m_popInitialMode = static_cast<PopInitMethod>(mode);
	addProTag(MKP);
	setOptType(MIN_OPT);

	v[param_sampleFre] = m_numDim * 2;

	Solution<CodeVInt>::allocateMemoryWB(m_numDim, m_numObj);
}

MultidimensionalKnapsack::~MultidimensionalKnapsack()
{
}

void MultidimensionalKnapsack::readProblem()
{
	size_t i;
	string Line;
	ostringstream oss;
	ifstream infile;
	oss << Global::g_arg[param_workingDir] << "Problem/Combination/MKP/data/" << m_fileName;
	infile.open(oss.str().c_str());
	if (!infile) {
		throw myException("read Multidimensional Knapsack data error");
	}
	infile >> Line;
	m_numDim = atoi(Line.c_str());
	infile >> Line;
	m_m = atoi(Line.c_str());
	infile >> Line;
	float temp = atof(Line.c_str());
	if (temp != 0) {
		vector<vector<double> > obj;
		obj.push_back(vector<double>(1, -temp));
		m_globalOpt.setGloObj(obj);
	}
	mv_p.resize(m_numDim);
	mv_b.resize(m_m);
	mvv_r.resize(m_m);
	for (i = 0; i < m_m; i++) {
		mvv_r[i].resize(m_numDim);
	}
	for (i = 0; i < m_numDim; i++)
	{
		infile >> mv_p[i];
		if (i == 0) {
			m_maxP = mv_p[i];
		}
		else if (m_maxP < mv_p[i]) {
			m_maxP = mv_p[i];
		}
	}
	for (i = 0; i < m_m; i++) {
		for (int j = 0; j < m_numDim; j++) {
			infile >> mvv_r[i][j];
		}
	}
	for (i = 0; i < m_m; i++) {
		infile >> mv_b[i];
	}
	infile.close();
	infile.clear();
}

int MultidimensionalKnapsack::invalidConstrainNum(VirtualEncoding &s_) {
	CodeVInt &s = dynamic_cast< CodeVInt&>(s_);
	int n = 0;
	double sum;
	for (int i = 0; i < m_m; i++) {
		sum = 0;
		for (int j = 0; j < m_numDim; j++) {
			sum += mvv_r[i][j] * s.m_x[j];
		}
		if (sum > mv_b[i]) n++;
	}
	return n;
}

ReturnFlag MultidimensionalKnapsack::evaluate_(VirtualEncoding &s_, bool rFlag, ProgramMode mode, bool flag2)
{
	CodeVInt &s = dynamic_cast< CodeVInt&>(s_);
	int m = invalidConstrainNum(s_);
	for (int i = 0; i<m_numObj; i++)
		s.m_obj[i] = 0;
	for (int n = 0; n<m_numObj; n++)
	{
		for (size_t i = 0; i<m_numDim; i++)
			s.m_obj[n] += mv_p[i] * s.m_x[i];
		s.m_obj[n] = s.m_obj[n] - m*m_maxP;
		s.m_obj[n] = -s.m_obj[n];
	}
	if (flag2) {
		if (rFlag)	m_evals++;
		if (mode == Program_Algorithm&&Global::msp_global->mp_problem) m_globalOpt.isFound(s.m_obj);

		if (Global::msp_global->mp_algorithm.get() != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()) return Return_Terminate;
		return Return_Normal;
	}
	return Return_Normal;
}

bool MultidimensionalKnapsack::getObjGlobalOpt(vector<double> &value) {
	if (m_globalOpt.flagGloObj()) {
		value = m_globalOpt[0].obj();
		return true;
	}
	return false;

}

bool MultidimensionalKnapsack::getObjGlobalOpt(vector<vector<double>> &value) {
	if (m_globalOpt.flagGloObj()) {
		value.clear();
		for (unsigned i = 0; i<m_globalOpt.getNumOpt(); i++)	value.push_back(m_globalOpt[i].obj());
		return true;
	}
	return false;

}

const Optima<CodeVInt> & MultidimensionalKnapsack::getGOpt()const {
	return m_globalOpt;
}

Optima<CodeVInt> & MultidimensionalKnapsack::getGOpt() {
	return m_globalOpt;
}

MultidimensionalKnapsack *MultidimensionalKnapsack::getTypePtr() {
	return this;
}
MultidimensionalKnapsack &MultidimensionalKnapsack::getTypeRef() {
	return *this;
}
bool MultidimensionalKnapsack::isGlobalOptKnown() {
	return m_globalOpt.flagGloObj();
}

void MultidimensionalKnapsack::setObjSet() {
	m_os.clear();
	if (!m_globalOpt.flagGloObj()) return;
	int num = m_globalOpt.getNumOpt();
	for (int i = 0; i < num; ++i) {
		m_os.push_back(&m_globalOpt[i].data().m_obj);
	}
}

double MultidimensionalKnapsack::getDistance(const VirtualEncoding &ss1, const VirtualEncoding &ss2, DistanceMode mode) {
	const CodeVInt& s1 = dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2 = dynamic_cast<const CodeVInt&>(ss2);
	double dis = 0;
	for (int i = 0; i<m_numDim; i++) {
		if (s1[i] != s2[i]) dis += 1;
	}
	return dis;
}

bool MultidimensionalKnapsack::isSame(const VirtualEncoding &ss1, const VirtualEncoding &ss2)
{
	const CodeVInt& s1 = dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2 = dynamic_cast<const CodeVInt&>(ss2);
	for (int i = 0; i<m_numDim; i++)
		if (s1[i] != s2[i])
			return false;
	return true;
}


void MultidimensionalKnapsack::initializeSolution(VirtualEncoding &result_, const int idx, const int maxId)
{
	CodeVInt& result = dynamic_cast< CodeVInt&>(result_);
	for (int i = 0; i < m_numDim; i++) {
		result[i] = Global::msp_global->getRandInt(0, 2);
	}
	if (!isValid(result))
		throw myException("error in @MultidimensionalKnapsack::initializeSolution() in MultidimensionalKnapsack.cpp");
}

bool MultidimensionalKnapsack::isValid(const VirtualEncoding &s) {
	const CodeVInt& s1 = dynamic_cast<const CodeVInt&>(s);
	for (int i = 0; i < m_numDim; i++) {
		if (s1.m_x[i] != 0 && s1.m_x[i] != 1)
			return false;
	}
	return true;
}