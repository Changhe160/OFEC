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
#ifndef FMAX_GLOBAL2_H
#define FMAX_GLOBAL2_H

#include "BenchmarkFunction.h"


class FMAX_global2 : public BenchmarkFunction
{
    public:
		FMAX_global2(ParamMap &v);
        FMAX_global2(const int rId, const int rDim,  string& rName);
        virtual ~FMAX_global2();
    protected:

        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};


#endif // FMAX_GLOBAL2_H
