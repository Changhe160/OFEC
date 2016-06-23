#ifndef	CONTINUOUS_PROBLEM
#define CONTINUOUS_PROBLEM
#include "problem.h"
#include "../Global/boundary.h"
#include "optimum.h"

#define CAST_PROBLEM_CONT dynamic_cast<ContinuousProblem*>(Global::msp_global->mp_problem.get())

class ContinuousProblem: public virtual Problem{
protected:
	double m_disAccuracy;
	BoundaryCont m_searchRange;
	Optima<CodeVReal> m_globalOpt; 
protected:
	ContinuousProblem():Problem(){}
	void parameterSetting(Problem * rP);
	void allocateMemory(const int numDim);
	void resizeDim(int num);
	void resizeObj(int num);
	void setObjSet();
public:
	ContinuousProblem(const int rId, const int rDimNumber, string rName, int numObj);
	bool isValid(const VirtualEncoding  &s);
	void validate(VirtualEncoding &s,SolutionValidation *mode=0);
	void validate_(vector<double> &s, SolutionValidation *mode = 0);
	void initializeSolution(VirtualEncoding &result,const int idx=0,const int maxId=0);
	void initializeSolution(const VirtualEncoding &ref,VirtualEncoding &result,double range);
	void initializePartSolution(VirtualEncoding &result,int begin,int end);
	ContinuousProblem& operator=(const ContinuousProblem &othr);
	virtual ~ContinuousProblem()=0;
	void getSearchRange(double &l,double&u, int i);
	double getDomainSize();
	bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2);
	void setDisAccuracy(double dis);
	inline double getDisAccuracy();

	bool getBoundaryFlag(int i);
	void setSearchRange(double l,double u);
	void setSearchRange(const vector<double> &l, const vector<double> &u);
	BoundaryCont & getSearchRange(){ return m_searchRange;}
	double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode);
	bool getObjGlobalOpt(vector<double> &opt);
	bool getObjGlobalOpt(vector<vector<double>> &opt);
	const Optima<CodeVReal> & getGOpt()const; 
	inline Optima<CodeVReal> & getGOpt();
	bool isGlobalOptKnown();
	const vector<pair<double, double>>& getObjRange();
	bool isParetoSet(const CodeVReal &s);  //check if s is dominated by the POS
	bool isGlobalOptFound();
	
};

inline double ContinuousProblem::getDisAccuracy(){
	return m_disAccuracy;
}

inline const Optima<CodeVReal> & ContinuousProblem::getGOpt()const{
	return m_globalOpt;
}

inline Optima<CodeVReal> & ContinuousProblem::getGOpt(){
	return m_globalOpt;
}
inline void ContinuousProblem::getSearchRange(double &l, double&u, int i){
	m_searchRange.getSearchRange(l, u, i);
}
inline bool  ContinuousProblem::isParetoSet(const CodeVReal &s){
	int num = m_globalOpt.getNumOpt();
	for (int i = 0; i < num; ++i){
		if (Compare_Better == m_globalOpt[i].compare_(s)) return false;
	}
	return true;
}

inline bool ContinuousProblem::isGlobalOptFound() {
	if (isGlobalOptKnown()) {
		if (m_globalOpt.isAllFound()) return true;
	}
	return false;
}
#endif