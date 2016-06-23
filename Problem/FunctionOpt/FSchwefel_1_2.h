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
#ifndef FSCHWEFEL_1_2_H
#define FSCHWEFEL_1_2_H

#include "BenchmarkFunction.h"


class FSchwefel_1_2 : public BenchmarkFunction
{
    public:
		FSchwefel_1_2(ParamMap &v);
        FSchwefel_1_2(const int rId, const int rDimNumber, string& rName);
        ~FSchwefel_1_2();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FSCHWEFEL_1_2_H
