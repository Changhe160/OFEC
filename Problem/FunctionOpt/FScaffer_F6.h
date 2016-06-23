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
// Created: 21 July 2011
// Last modified:
#ifndef FSCAFFER_F6_H
#define FSCAFFER_F6_H

#include "BenchmarkFunction.h"


class FScaffer_F6 : public BenchmarkFunction
{
    public:
		FScaffer_F6(ParamMap &v);
        FScaffer_F6(const int rId, const int rDim, string& rName);
        ~FScaffer_F6();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FSCAFFER_F6_H
