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
// Created: 20 July 2011
// Last modified:
#ifndef BENCHMARKFUNCTION_H
#define BENCHMARKFUNCTION_H

#include "../ContinuousProblem.h"
#include "../../Utility/matrix.h"
class BenchmarkFunction : public ContinuousProblem
{

public:
    BenchmarkFunction(const int rId, const int rDimNumber,  string rName, int numObj=1);       
    ~BenchmarkFunction();
    void setTranslation(const double *rTrans);
    void setRotation(const double * const * rRot);
    void setBias(double rBias);
    void setScale(double rScale);
    void setRotationFlag(bool rFlag);
    void setTranlationFlag(bool rFlag);
    
    double * getTranslation();
    Matrix * getRotation();
    double getConditionNumber();
    void setConditionNumber(double rC);
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true);
	int & idx(){ return m_idx; }
	int & gidx(){ return m_gidx; }
	int &tidx(){ return m_tidx; };
protected:
    double *mp_translation;         // shifted optimum positions for f(x)
    bool m_scaleFlag,m_rotationFlag,m_translationFlag,m_biasFlag;
    double m_scale,m_bias;
    double m_conditionNumber;       // generation ratation matrix for functions with rotation property
    Matrix *mp_rotationMatrix;      // rotation matrix
	int m_gidx = 0, m_idx = 0,m_tidx=0;	//index of a group it belongs to, index in a group, and index of a tree it is 

    virtual void freeMemory();
    virtual void allocateMemory(const int rDim);
    BenchmarkFunction& operator=(const BenchmarkFunction & rBF);

    void zeroTranslation();
    bool isZeroTranslation();
    virtual void transform(double * x);

    virtual void initialize(){};
    virtual bool loadTranslation();
    virtual bool loadRotation();
    void setGlobalOpt(int i=0,vector<double>* opt=0, double *tran=0);
    void setOriginalGlobalOpt(int i=0,vector<double>* opt=0);	
	virtual void evaluate__(double const *x,vector<double>& obj)=0;
public:
	Optima<CodeVReal,Solution<CodeVReal>> m_originalGlobalOpt;
};

#endif // BENCHMARKFUNCTION_H
