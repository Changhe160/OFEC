//#pragma once
#ifndef TERMINATION_H
#define TERMINATION_H
#include <vector>
#include "../Utility/TypeVar/typeVar.h"
#include "../Utility/TypeList/Typelist.h"


using namespace std;

//terminate when the global opt. is found
class Termination 
{
public:
	Termination(ParamMap &v) {}
	virtual ~Termination() {}
	
	virtual bool ifTerminating();
	bool ifTerminated() { return m_isTerminated; }
	void setTermTrue() { m_isTerminated = true; }
	
protected:
	bool m_isTerminated=false;
};

//terminate when the maximum number of iterations is reached
class TermMaxGen :public Termination {
protected:
	int m_maxIter=1000;
public:
	using Termination::ifTerminating;
	TermMaxGen(ParamMap &v) :Termination(v) {
		if (v.find(param_maxIter) != v.end()) m_maxIter = v[param_maxIter];
	}
	bool ifTerminating(int value) {
		// Assume that vp[0] is a main population that records the number of iterations since the run starts
		
#if defined OFEC_DEMON
		return Termination::ifTerminating();
#else
		if (Termination::ifTerminating()) return true;
#endif
	
		if (value >= m_maxIter) return true;
		
		return false;
	}
};

//terminate when the maximum number of evaluations is reached
class TermMaxFes :public Termination {
protected:
	int m_maxFes=10000;
public:
	TermMaxFes(ParamMap &v) :Termination(v){
	
		if (v.find(param_maxEvals) != v.end()) m_maxFes = v[param_maxEvals];
	
	}

	bool ifTerminating();
};

//terminate when the best solution remains over a number of successive iterations
class TermBest :public Termination {
protected:
	int m_maxSucIter=200;
	int m_sucIter = 0;
	double m_preBest,m_curBest;
public:
	using Termination::ifTerminating;
	TermBest(ParamMap &v) :Termination(v){
		if (v.find(param_maxSucIter) != v.end()) m_maxSucIter = v[param_maxSucIter];
	}
	
	void initialize(const vector<TypeVar> *value) {	
		m_preBest = m_curBest = value->at(0);		
	}
	bool ifTerminating(vector<double> &value) {
#if defined OFEC_DEMON
		return Termination::ifTerminating();
#else
		if (Termination::ifTerminating()) return true;
#endif

		if (m_sucIter > m_maxSucIter) return true;

		return false;
	}
	void countSucIter(vector<double> &value) {
		m_curBest = value[0];

		if (m_curBest != m_preBest) {
			m_preBest = m_curBest;
			m_sucIter = 0;
		}
		else {
			m_sucIter++;
		}
	}
};

//terminate when the average objectives changes less than a threshold value over a number of successive iterations
class TermMean :public Termination {
protected:
	int m_maxSucIter = 200;
	int m_sucIter = 0;
	double m_preMean, m_curMean,m_epsilon=1.E-2;
public:
	using Termination::ifTerminating;
	TermMean(ParamMap &v) :Termination(v) {
		if (v.find(param_epsilon) != v.end()) m_epsilon = v[param_epsilon];	
		if (v.find(param_maxSucIter) != v.end()) m_maxSucIter = v[param_maxSucIter];
	}
	void initialize(double value) {
		m_preMean = m_curMean = value;
	}
	bool ifTerminating(double value);
	void countSucIter(double value);
};

//terminate when the variance of objective ls less than a small value
class TermVar :public Termination {
protected:
	double m_epsilon = 1E-3;
public:
	using Termination::ifTerminating;
	TermVar(ParamMap &v) :Termination(v) {
		if (v.find(param_epsilon) != v.end()) m_epsilon = v[param_epsilon];
}

	bool ifTerminating(double mean, double var) {
#if defined OFEC_DEMON
		return Termination::ifTerminating();
#else
		if (Termination::ifTerminating()) return true;
#endif
	
		if (mean < 1) m_epsilon *= m_epsilon;
		if (var < m_epsilon) return true;
		return false;
	}
	

};

typedef LOKI_TYPELIST_6(Termination, TermMaxFes, TermMaxGen, TermBest, TermMean, TermVar) TermList;


#endif // !TERMINATION_H

