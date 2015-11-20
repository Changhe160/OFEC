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

/************************************************************************
parameter estimation for frequency-modulated (FM) sound waves is to estimate
the parameters of a FM synthesizer.

S.~Das and P.~N. Suganthan, ``Problem definitions and evaluation criteria for
  cec 2011 competition on testing evolutionary algorithms on real world
  optimization problem,'' Dept. of Electronics and Telecommunication Engg.,
  Jadavpur University, Kolkata, India, Tech. Rep., 2011.

*****************************************************************************************/



// Created: 21 July 2011
// Last modified:
#ifndef FPAREST_FMSOUNDWAVES_H
#define FPAREST_FMSOUNDWAVES_H


#include "../FunctionOpt/BenchmarkFunction.h"

class FParEst_FMSoundWaves : public BenchmarkFunction
{
    public:
        FParEst_FMSoundWaves(ParamMap &v);
        FParEst_FMSoundWaves(const int rId,  const int rDim, string& rName);
        virtual ~FParEst_FMSoundWaves();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FPAREST_FMSOUNDWAVES_H
