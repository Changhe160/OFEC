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
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#ifndef FBEASLEY_H
#define FBEASLEY_H

#include "BenchmarkFunction.h"

class FBeasley : public BenchmarkFunction
{ //Beasley's F4 function
    public:
		FBeasley(ParamMap &v);
        FBeasley(const int rId, const int rDimNumber, string& rName);
        virtual ~FBeasley();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FBEASLEY_H