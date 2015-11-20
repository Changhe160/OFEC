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
#ifndef FMAX_GLOBAL5_H
#define FMAX_GLOBAL5_H

#include "BenchmarkFunction.h"


class FMAX_global5 : public BenchmarkFunction
{
    public:
		FMAX_global5(ParamMap &v);
        FMAX_global5(const int rId,  const int rDim, string& rName);
        virtual ~FMAX_global5();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FMAX_GLOBAL5_H
