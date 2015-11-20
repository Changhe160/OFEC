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
#ifndef FACKLEY_H
#define FACKLEY_H

#include "BenchmarkFunction.h"


class FAckley : public BenchmarkFunction
{
    public:
		FAckley(ParamMap &v);
        FAckley(const int rId, const int rDimNumber, string& rName);
        virtual ~FAckley();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FACKLEY_H
