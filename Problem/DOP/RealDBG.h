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

#ifndef REALDBG_H
#define REALDBG_H

#include "DynamicContinuous.h"
#include "../../Utility/matrix.h"

class RealDBG : public DynamicContinuous
{
    protected:
        bool m_prediction;						// the next change of function can be predictable or not
        Matrix *mp_rotationMatrix;				// orthogonal rotation matrixes for each function
        int ***mppp_rotationPlanes;					// save the planes rotated during one periodicity

    public:
        RealDBG(const int rId, const int rDimNumber, const int rNumPeaks, const int numObj=1);
        virtual ~RealDBG()=0;

        ///TODO: not 100% sure whether it should be derived by CompositionDBG and RotationDBG
        virtual bool setPeriod(const int rPeriod);

        RealDBG & operator=(const RealDBG & rP);
        void reset();
		void reinitialize();
    protected:
        void correctSolution(double *x);
        void heightStandardChange();
        void positionStandardChange(double angle);

        virtual void randomChange(){};
        virtual void smallStepChange(){};
        virtual void largeStepChange(){};
        virtual void recurrentChange(){};
        virtual void chaoticChange(){};
        virtual void recurrentNoisyChange(){};


        void parameterSetting(Problem * rP);
        virtual void  freeMemory();
        double  standardChange(const ChangeType T, const double min, const double max);
        virtual void allocateMemory(const int rDimNum, const int rPeaks);

        void initialize(const ChangeType rT, const bool rFlagDimChange, const bool rFlagNumPeakChange,
			const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage);
        void restoreInfor();
        virtual void changeDimension(){};
        virtual void changeNumPeaks(){};

};

#endif // REALDBG_H
