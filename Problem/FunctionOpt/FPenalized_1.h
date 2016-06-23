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
#ifndef FPENALIZED_1_H
#define FPENALIZED_1_H

#include "BenchmarkFunction.h"


class FPenalized_1 : public BenchmarkFunction
{
    public:
		FPenalized_1(ParamMap &v);
        FPenalized_1(const int rId, const int rDimNumber, string& rName);
        virtual ~FPenalized_1();
    protected:
        void initialize();
        double u(double x, double a, double k, double m)const;
        void evaluate__(double const *x,vector<double>& obj);
    private:
};
#endif // FPENALIZED_1_H
