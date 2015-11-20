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
/*
J. Branke, Memory enhanced evolutionary algorithms for changing optimization
problems, Proc. of the 1999 Congr. on Evol. Comput, pp. 1875-1882, 1999.
*/
// Created: 16 May 2013
// Last modified:
#ifndef MOVINGPEAK_H
#define MOVINGPEAK_H

#include "DynamicContinuous.h"

class MovingPeak : public DynamicContinuous
{
private:
    int m_F;
    double m_vlength ; /* distance by which the peaks are moved, severity */

    /* lambda determines whether there is a direction of the movement, or whether
    they are totally random. For lambda = 1.0 each move has the same direction,
    while for lambda = 0.0, each move has a random direction */
    double m_lambda;
    int m_useBasisFunction; /* if set to 1, a static  landscape (basis_function) is included in the fitness evaluation */
    int m_calculateRightPeak ; /* saves computation time if not needed and set to 0 */
    double  m_standardHeight;
    /* width chosen randomly when standardwidth = 0.0 */
    double  m_standardWidth;

    double *mp_shift;
    int *mp_coveredPeaks;    /* which peaks are covered by the population ? */
    double **mpp_prevMovement;/* to store every peak's previous movement */
	

private:

	bool readData();

    /* the following basis functions are provided :*/
    double constantBasisFunc(const double *gen);
    double fivePeakBasisFunc(const double *gen);
    /* the following peak functions are provided: */
    double peakFunction1(const double *gen, int peak_number);
    double peakFunctionCone (const double *gen, const int &peak_number);
    double peakFunctionHilly (const double *gen, int peak_number);
    double peakFunctionTwin (const double  *gen, int peak_number);
    double functionSelection(const double  *gen, const int &peak_number);
    double dummyEval (const double *gen);
    void initialize();
	

protected:
    void currentPeakCalc (const double *gen);
    virtual void  freeMemory();
    virtual void allocateMemory(const int rDimNum, const int rPeaks);
    virtual void randomChange();
    ///TODO the flowing change types are not implemented
    virtual void smallStepChange(){randomChange();};
    virtual void largeStepChange(){randomChange();};
    virtual void recurrentChange(){randomChange();};
    virtual void chaoticChange(){randomChange();};
    virtual void recurrentNoisyChange(){randomChange();};
    virtual void changeDimension(){randomChange();};
    virtual void changeNumPeaks();
	void parameterSetting(Problem * rP);
		//16/05/2013
	virtual void updatePeakQaulity();
	virtual void calculateAssociateRadius();
	void setSeverity();
	
public:
    MovingPeak(const int rId, const int rDimNumber, const int rNumPeaks, double const rChangingRatio,const bool rFlagDimChange=false,\
		const bool rFlagNumPeakChange=false,const int peakNumChangeMode=1,const bool flagNoise=false, const bool flagTimelinkage=false);
    MovingPeak(ParamMap &v);
    ~MovingPeak();
	
    void changeStepsizeRandom () ;
    void changeStepsizeLinear();
    int getRightPeak();
    void setVlength(const double s);
    void reset();
	MovingPeak &operator=(MovingPeak &other);
	//void reinitialize();
	double getVLength();
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true);
};

#endif // MOVINGPEAK_H
