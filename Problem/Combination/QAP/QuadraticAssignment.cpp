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


#include "QuadraticAssignment.h"
#include "../../../Global/global.h"
#include<string>

using namespace std;

QuadraticAssignment::QuadraticAssignment(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1), m_globalOpt(m_numDim, m_numObj)
{
	v[param_numObj] = 1;
	mvv_flow.resize(m_numDim);
	mvv_distance.resize(m_numDim);
	for (int i = 0; i < m_numDim; i++) {
		mvv_flow[i].resize(m_numDim);
		mvv_distance[i].resize(m_numDim);
	}
	m_fileName = (string)v[param_dataFile1];
	readProblem();

	m_searchRange.setBoundary(0, (v[param_numDim]) - 1);
	int mode = int(v[param_populationInitialMethod]);
	m_popInitialMode = static_cast<PopInitMethod>(mode);
	addProTag(QAP);
	setOptType(MIN_OPT);

	v[param_sampleFre] = m_numDim * 2;

	Solution<CodeVInt>::allocateMemoryWB(m_numDim, m_numObj);
}

QuadraticAssignment::~QuadraticAssignment()
{
}

bool QuadraticAssignment::isValid(const VirtualEncoding &s1)
{
	const CodeVInt& s = dynamic_cast<const CodeVInt&>(s1);

	for (int i = 0; i<m_numDim; i++) {  //judge the range		
		if ((s[i])<m_searchRange.m_lower || (s[i])>m_searchRange.m_upper) return false;
	}
	vector<int> flag(m_numDim, 0);  //judge whether has the same gene
	int temp;
	for (int i = 0; i<m_numDim; i++)
	{
		temp = s[i];
		flag[temp] = 1;
	}
	for (int i = 0; i<m_numDim; i++)
		if (flag[i] == 0)
			return false;
	return true;
}

void QuadraticAssignment::readProblem()
{
	size_t i;
	string Line;
	char *Keyword = 0;
	const char *Delimiters = " ():=\n\t\r\f\v\xef\xbb\xbf";
	ostringstream oss;
	ifstream infile;
	oss << Global::g_arg[param_workingDir] << "Problem/Combination/QAP/data/" << m_fileName;
	infile.open(oss.str().c_str());
	if (!infile) {
		throw myException("read Quadratic Assignment data error");
	}
	char *savePtr;
	while (getline(infile, Line))
	{
		if (!(Keyword = gStrtok_r((char*)Line.c_str(), Delimiters, &savePtr)))
			continue;
		for (i = 0; i<strlen(Keyword); i++)
			Keyword[i] = toupper(Keyword[i]);
		if (!strcmp(Keyword, "DIM"))
		{
			char *token = gStrtok_r(0, Delimiters, &savePtr);
			m_numDim = atoi(token);
		}
		else if (!strcmp(Keyword, "FLOW"))
		{
			for (int n = 0; n < m_numDim; n++)
				for (i = 0; i < m_numDim; i++)
					infile >> mvv_flow[n][i];
		}
		else if (!strcmp(Keyword, "DISTANCE"))
		{
			for (int n = 0; n < m_numDim; n++)
				for (i = 0; i < m_numDim; i++)
					infile >> mvv_distance[n][i];
		}
		else if (!strcmp(Keyword, "OPT_OBJ"))
		{
			char *token = gStrtok_r(0, Delimiters, &savePtr);
			vector<vector<double> > obj;
			obj.push_back(vector<double>(1, atof(token)));
			m_globalOpt.setGloObj(obj);
		}
		else if (!strcmp(Keyword, "OPT_SOLUTION"))
		{
			int temp, n=0;
			infile >> Line;
			for (i = 1; n < m_numDim; i += 2) {
				temp = Line[i] - '0';
				--temp;
				m_globalOpt[0].data().m_x[n++] = TypeVar(temp);
			}
			m_globalOpt.setFlagLocTrue();
		}
	}
	infile.close();
	infile.clear();
}


ReturnFlag QuadraticAssignment::evaluate_(VirtualEncoding &s_, bool rFlag, ProgramMode mode, bool flag2)
{
	CodeVInt &s = dynamic_cast< CodeVInt&>(s_);

	for (int i = 0; i<m_numObj; i++)
		s.m_obj[i] = 0;
	int row, col;
	for (int n = 0; n<m_numObj; n++)
	{
		for (size_t i = 0; i<m_numDim; i++)
		{
			row = s.m_x[i];
			for (size_t j = 0; j < m_numDim; j++) {
				col = s.m_x[j];
				s.m_obj[n] += mvv_distance[i][j] * mvv_flow[row][col];
			}
		}
	}
	if (flag2) {
		if (rFlag)	m_evals++;
		if (mode == Program_Algorithm&&Global::msp_global->mp_problem) m_globalOpt.isFound(s.m_obj);

		if (Global::msp_global->mp_algorithm.get() != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()) return Return_Terminate;
		return Return_Normal;
	}
	return Return_Normal;
}

bool QuadraticAssignment::getObjGlobalOpt(vector<double> &value) {
	if (m_globalOpt.flagGloObj()) {
		value = m_globalOpt[0].obj();
		return true;
	}
	return false;

}

bool QuadraticAssignment::getObjGlobalOpt(vector<vector<double>> &value) {
	if (m_globalOpt.flagGloObj()) {
		value.clear();
		for (unsigned i = 0; i<m_globalOpt.getNumOpt(); i++)	value.push_back(m_globalOpt[i].obj());
		return true;
	}
	return false;

}

const Optima<CodeVInt> & QuadraticAssignment::getGOpt()const {
	return m_globalOpt;
}

Optima<CodeVInt> & QuadraticAssignment::getGOpt() {
	return m_globalOpt;
}

QuadraticAssignment *QuadraticAssignment::getTypePtr() {
	return this;
}
QuadraticAssignment &QuadraticAssignment::getTypeRef() {
	return *this;
}
bool QuadraticAssignment::isGlobalOptKnown() {
	return m_globalOpt.flagGloObj();
}

void QuadraticAssignment::setObjSet() {
	m_os.clear();
	if (!m_globalOpt.flagGloObj()) return;
	int num = m_globalOpt.getNumOpt();
	for (int i = 0; i < num; ++i) {
		m_os.push_back(&m_globalOpt[i].data().m_obj);
	}
}

double QuadraticAssignment::getDistance(const VirtualEncoding &ss1, const VirtualEncoding &ss2, DistanceMode mode) {
	const CodeVInt& s1 = dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2 = dynamic_cast<const CodeVInt&>(ss2);
	double dis = 0;
	for (int i = 0; i<m_numDim; i++) {
		if (s1[i] != s2[i]) dis += 1;
	}
	return dis;
}

bool QuadraticAssignment::isSame(const VirtualEncoding &ss1, const VirtualEncoding &ss2)
{
	const CodeVInt& s1 = dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2 = dynamic_cast<const CodeVInt&>(ss2);
	for (int i = 0; i<m_numDim; i++)
		if (s1[i] != s2[i])
			return false;
	return true;
}


void QuadraticAssignment::initializeSolution(VirtualEncoding &result_, const int idx, const int maxId)
{
	CodeVInt& result = dynamic_cast< CodeVInt&>(result_);
	vector<int> temp;
	int i, pos, num = result.m_x.size();
	for (i = 0; i<num; i++)
		temp.push_back(int(i));
	for (i = 0; i<num; i++)
	{
		pos = int((num - 1 - i)*Global::msp_global->mp_uniformAlg->Next());
		result[i] = temp[pos];
		temp[pos] = temp[num - 1 - i];
	}
	if (!isValid(result))
		throw myException("error in @QuadraticAssignment::initializeSolution() in QuadraticAssignment.cpp");
}