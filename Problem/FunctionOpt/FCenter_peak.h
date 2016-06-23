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
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#ifndef FCENTER_PEAK_H
#define FCENTER_PEAK_H

#include "BenchmarkFunction.h"
//Ursem F4 in R. K. Ursem, ¡°Multinational evolutionary algorithms,¡± in Proc. Congr.
//Evol. Comput. (CEC), vol. 3. 1999, pp. 1633¨C1640.

class FCenter_peak : public BenchmarkFunction
{//Ursem F4
    public:
		FCenter_peak(ParamMap &v);
        FCenter_peak(const int rId, const int rDimNumber, string& rName);
        virtual ~FCenter_peak();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FCENTER_PEAK_H