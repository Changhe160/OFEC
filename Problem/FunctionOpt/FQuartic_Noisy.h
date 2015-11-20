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
#ifndef FQUARTIC_NOISY_H
#define FQUARTIC_NOISY_H

#include "BenchmarkFunction.h"

class FQuartic_Noisy : public BenchmarkFunction
{
    public:
		FQuartic_Noisy(ParamMap &v);
        FQuartic_Noisy(const int rId, const int rDimNumber, string& rName);
        virtual ~FQuartic_Noisy();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FQUARTIC_NOISY_H
