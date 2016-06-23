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
#ifndef FSIX_HUMP_H
#define FSIX_HUMP_H

#include "BenchmarkFunction.h"
//Stoean, C.; Preuss, M.; Stoean, R.; Dumitrescu, D., "Multimodal Optimization by Means of a Topological Species Conservation Algorithm," Evolutionary Computation, IEEE Transactions on , vol.14, no.6, pp.842,864, Dec. 2010
//doi: 10.1109/TEVC.2010.204166

class FSix_humpCamelBack : public BenchmarkFunction
{
    public:
		FSix_humpCamelBack(ParamMap &v);
        FSix_humpCamelBack(const int rId, const int rDimNumber, string& rName);
        virtual ~FSix_humpCamelBack();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FSIX_HUMP_H