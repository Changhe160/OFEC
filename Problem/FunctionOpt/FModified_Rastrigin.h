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
#ifndef FMODIFIED_RASTRIGIN
#define FMODIFIED_RASTRIGIN

#include "BenchmarkFunction.h"

class FModified_Rastrigin : public BenchmarkFunction
{ 
    public:
		FModified_Rastrigin(ParamMap &v);
        FModified_Rastrigin(const int rId, const int rDimNumber, string& rName);
        virtual ~FModified_Rastrigin();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
		vector<double> m_k;
};

#endif 