//#pragma once
#ifndef OneMax_H
#define OneMax_H
#include "../../problem.h"
#include "../../optimum.h"
#define CAST_OneMax dynamic_cast<OneMax *>(Global::msp_global->mp_problem.get())
class OneMax:public Problem {
protected:
	Optima<CodeVInt> m_globalOpt;
	virtual void setObjSet() = 0;
public:
	OneMax(ParamMap& v);
	OneMax(const int rId, const int rDimNumber, string rName, const int numObj = 1);
	~OneMax() {};
	
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag = true) = 0;

	bool isValid(const VirtualEncoding &s) = 0;
	void validate(VirtualEncoding &s, SolutionValidation *mode = 0) {};
	void initializeSolution(VirtualEncoding &result, const int idx = 0, const int maxId = 0);
	void initializeSolution(const VirtualEncoding &ref, VirtualEncoding &result, double range) {};
	void initializePartSolution(VirtualEncoding &result, int begin, int end) {};

	double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode);
	bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2);
	bool isGlobalOptKnown() {
		return true;
	};

	bool getObjGlobalOpt(vector<double> &opt) {
		opt = m_globalOpt[0].data().m_obj;
		return true;
	};
	bool getObjGlobalOpt(vector<vector<double>> &value) {
		
		value.clear();
		for (unsigned i = 0; i<m_globalOpt.getNumOpt(); i++)	value.push_back(m_globalOpt[i].obj());
		return true;
		
	}
	const Optima<CodeVInt> & getGOpt()const {
		return m_globalOpt;
	}

	Optima<CodeVInt> & getGOpt() {
		return m_globalOpt;
	}
	
};
#endif // !OneMax_H

