/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
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
#ifndef ROTATIONDBG_H
#define ROTATIONDBG_H

#include "RealDBG.h"


class RotationDBG : public RealDBG
{
public:
	RotationDBG(const int rId, const int rDimNumber,  const int rNumPeaks, \
             const ChangeType rT,double const rChangingRatio, const bool rFlagDimChange,const bool rFlagNumPeakChange,\
			const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage );
	RotationDBG(ParamMap &v);
	virtual ~RotationDBG();    
	RotationDBG& operator=(const RotationDBG &);
	virtual void  setWidth(const double w);
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true);
protected:
	void parameterSetting(Problem * rP);
    void widthStandardChange();
	virtual void randomChange();
    virtual void smallStepChange();
    virtual void largeStepChange();
    virtual void recurrentChange();
    virtual void chaoticChange();
    virtual void recurrentNoisyChange();
    virtual void changeDimension();
    virtual void changeNumPeaks();
	void initialize(const ChangeType T,double const rChangingRatio, const bool rFlagDimChange=false,const bool rFlagNumPeakChange=false,
	const int peakNumChangeMode=1,const bool flagNoise=false, const bool flagTimelinkage=false);
	
};

#endif // ROTATIONDBG_H
