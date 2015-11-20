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
// Created: 11 May 2011
// Last modified:

/// class for composition dynamic benchmark generator
///C. Li, S. Yang, T. T. Nguyen, E. L. Yu, X. Yao, Y. Jin, H.-G. Beyer, and
///P. N. Suganthan, ''Benchmark Generator for CEC'2009 Competition on
///Dynamic Optimization,'' Technical Report 2008, Department of Computer
///Science, University of Leicester, U.K., 2008.

#ifndef COMPOSITIONDBG_H
#define COMPOSITIONDBG_H

#include"RealDBG.h"

class CompositionDBG : public RealDBG
{
	enum ProblemTag{Sphere=0,Rastrigin,Griewank,Ackley,Weierstrass};
private:
	vector<Boundary> mv_comBoundary;				// boundary of component functions
	double * mp_convergeSeverity;						// severity of converge range for each function
	double * mp_stretchSeverity;							// severity of stretching original function, greater than 1 for stretch
														// less than 1 for compress the original function
	double m_heightNormalizeSeverity;					// constant number for noralizing all basic function with similar height
	ProblemTag *mp_comFunction;						// which basic function used to compose the composition function

    double *mp_fmax;
	static const int msc_numComFuns=5;					// number of basic functions
	static boost::thread_specific_ptr<ComDBGFuncID> ms_funID;
private:
	
	void setComBoundary();
	void getComBoundary(const ProblemTag &f,double &l,double &u,const int rDimIdx=0)const;
	void initialize(const ChangeType T, const ComDBGFuncID rF,  const double rChangingRatio,const bool rFlagDimChange=false,\
		const bool rFlagNumPeakChange=false,const int peakNumChangeMode=1,const bool flagNoise=false, const bool flagTimelinkage=false);

public:
    static const int msc_numFuns=5;                     // number of functions in comDBG system
	CompositionDBG(const int rId, const int rDimNumber,  const int rNumPeaks,const ChangeType rT, const ComDBGFuncID rF, \
		const double rChangingRatio,const bool rFlagDimChange,const bool rFlagNumPeakChange,\
				const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage);
	CompositionDBG(ParamMap &v);

    virtual ~CompositionDBG();

	CompositionDBG &operator =(const CompositionDBG &);
	void setRotationMatrix();					//randomly generate rotation matrx for each basic function
	void setCovergeSevrity(const double* cs);
	void setStretchSeverity();
	void setBasicFunction(const ProblemTag *bf);		//component functions to compose the search space	
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true);
protected:
	void parameterSetting(Problem * rP);
	virtual void  freeMemory();

	virtual void randomChange();
    virtual void smallStepChange();
    virtual void largeStepChange();
    virtual void recurrentChange();
    virtual void chaoticChange();
    virtual void recurrentNoisyChange();
    virtual void allocateMemory(const int rDimNum, const int rPeaks);
    virtual void changeDimension();
    virtual void changeNumPeaks();
	
private:
	void correctSolution(const ProblemTag &f, double *x);					// make x within search range after rotation											// basic five functions
	double fSphere(double *x);
	double fRastrigin(double *x);
	double fWeierstrass(double *x);
	double fGriewank(double *x);
	double fAckley(double *x);
	double selectFun(const ProblemTag &f,double *x);

};
#endif // COMPOSITIONDBG_H





