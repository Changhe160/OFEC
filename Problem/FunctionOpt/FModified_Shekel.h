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
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#ifndef FMODIFIED_SHEKEL_H
#define FMODIFIED_SHEKEL_H

#include "BenchmarkFunction.h"


class FModified_Shekel : public BenchmarkFunction
{
    public:
		FModified_Shekel(ParamMap &v);
        FModified_Shekel(const int rId, const int rDimNumber, string& rName);
        virtual ~FModified_Shekel();
    protected:
        void initialize();
         void evaluate__(double const *x,vector<double>& obj);
    private:
		double m_a[8][5];
		double m_c[8];
};

#endif // FMODIFIED_SHEKEL_H