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
#ifndef FROSENBROCK_H
#define FROSENBROCK_H

#include "BenchmarkFunction.h"


class FRosenbrock : public BenchmarkFunction
{
    public:
		FRosenbrock(ParamMap &v);
        FRosenbrock(const int rId, const int rDimNumber, string& rName);
        virtual ~FRosenbrock();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FROSENBROCK_H
