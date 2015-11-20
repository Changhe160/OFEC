/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
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
// Created: 21 July 2011
// Last modified:

#ifndef SINGLEOBJSTA_H
#define SINGLEOBJSTA_H

class Global;
#include "../Utility/TypeVar/typeVar.h"

class mSingleObj
{
    public:      
        static mSingleObj* getSingleObj();
        static void initialize(ParamMap &v);
      
        virtual ~mSingleObj();
        mSingleObj& operator=(const mSingleObj& other);
        double getConvergeSpeed2Best();
        double getConvergeSpeed2GOpt();
        double getConvergeSpeed();
        void getStatInfor(double &mean, double &var, double &worst, double &best);
        double getSucRate();
        double getPerformance();
		double getBestSoFar(int );
		const vector< vector<double>> & getGOpt() const { return mpp_gOpt; }
        void record(Global *,double rObj);
        void calculatePerformance();
        void addGOpt(int runId,double rGOptObj);
        void setFileName(stringstream &rName);
		void setFileName(ParamMap &v);
        void setProgrOutputFlag(bool rFlag);
        void outputResult();
        stringstream m_fileName,m_algPar,m_proPar;

        void setConvgProgMode(ConvgProgeMode rMode);
        void setAlgParameter(stringstream & rPar);
        void resetGOptIndx();
		void setAbsoluteErrFlag(bool flag=true);
		void setAccuracy(double acc);
		void setCompareType(Compare comp);
		void setProParameter(stringstream & rPar);
		static	void deleteSingleObj();
    protected:
        void calculateConvergenceTime();
		static unique_ptr<mSingleObj>  msp_perf;
        mSingleObj(ParamMap &v);
		virtual void calculateKeyParam();
		void calculateNumRecords();
    protected:
        vector< vector<double>> mpp_data;                                          //best obj value of each sampling record
        vector< vector<double>> mpp_gOpt;                                          //obj value of global optima
        vector<int> m_gOptIdx;                                              // index of current global optimum in mpp_gOpt
        double m_mean,m_var,m_largest,m_smallest;                   //the mean, variance(over changes), largest, smallest value in all runs
        double m_speed2Best,m_speed2gOpt,m_speed;                   // the convergence speed to the algorithm's best, the gOpt, and the fitness decrease/incscrease per evaluation
        double m_performance,m_perfVar;                             // the algorithm's performance defined in the cpp file and variance of performance over all runs 
        double m_sucRate;                                           // the success rate of achieving a certain accuracy level
        int m_numEvals2Suc;                                         // the number of evaluations needed to achieve the given succuss rate

        vector<double> m_bestSoFar;                                        // the best value since last record
        vector<vector<int>>mpp_convergeTime;                                    // the number of fitness evaluations when the algorithm converges
        int m_numRecords;                                          // the number of records needed to be recorded during one run
        int m_numChanges;                                            // the number of changes during one run
        int m_recordsPerChange;                                    // the number of records per change for DOPs (noted: one change in static optimization probelm)
        int m_numEvals2Converge;                                    // the number of evaluations till the algorithm converges

        bool m_progrOutputFlag;
        ConvgProgeMode m_convgMode;                                 // the mode of the convergence graph
		bool m_absoluteErr;		
		double m_offlineErrOverChanges, m_offlineErrOverChangesVar;	// offline error, average all fitness errors at each fitness evaluation(due to space issue with the method used, it is
		double m_offlineErrOverRuns,m_offlineErrOverRunsVar;		// evaluated every Global::g_arg[param_sampleFre] fes). The two ways differs only with the variance calculation(i.e., the same error value but different standard variance).
		double m_meanOverRuns,m_meanOverRunsVar;
		double m_accuracy;
		int m_changeFre;
		Compare m_comp;
		double m_avgEvals,m_avgCevals,m_avgTevals;
};

#endif // SINGLEOBJSTA_H
