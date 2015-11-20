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
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#ifndef FFIVE_HILLS_H
#define FFIVE_HILLS_H

#include "BenchmarkFunction.h"
//Ursem F3 in R. K. Ursem, ¡°Multinational evolutionary algorithms,¡± in Proc. Congr.
//Evol. Comput. (CEC), vol. 3. 1999, pp. 1633¨C1640.

class FFive_hills : public BenchmarkFunction
{
    public:
		FFive_hills(ParamMap &v);
        FFive_hills(const int rId, const int rDimNumber, string& rName);
        virtual ~FFive_hills();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FFIVE_HILLS_H