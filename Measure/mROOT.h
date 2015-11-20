#ifndef ROOT_H
#define ROOT_H


template<typename > class Solution;

#include "../Global/global.h"
#include "../Global/solution.h"

class mROOT{  // Robust Optimization Over Time, Only for DOPs with a single objective
private:
	unsigned m_timeWindow;							//the number of changes to observe
	vector<int> m_numChanges;							//the current number of changes
	vector<deque<vector< Solution<> >>> mq_solution;		//to store solutions within one time window
	vector<vector<double>> mv_root;						//to store root value of each run
	vector<double> mv_avgRoot;					//to store avg root values of all runs
	double m_mean,m_std;
private:
	static unique_ptr<mROOT> msp_root;
private:
	mROOT(ParamMap &v);
public:
	static mROOT * getROOT(){
		return msp_root.get();
	}
	~mROOT();
	static bool initilizeROOT(ParamMap &v);
	static void deleteRoot();
	bool updateObjValue(Global*);
	double getMeanRoot();
	double getROOTStd();
	void record(const Solution<> & chr,Global *glob);
	void calculateMean();
};

#endif