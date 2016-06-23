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
#ifndef FNONCONT_RASTRIGIN_H
#define FNONCONT_RASTRIGIN_H

#include "BenchmarkFunction.h"

class FNoncont_Rastrigin : public BenchmarkFunction
{
    public:
		FNoncont_Rastrigin(ParamMap &v);
        FNoncont_Rastrigin(const int rId, const int rDim, string& rName);
        virtual ~FNoncont_Rastrigin();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FNONCONT_RASTRIGIN_H
