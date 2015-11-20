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
*  Paper: Minimization of Keane¡¯s Bump Function by the Repulsive Particle Swarm and the 
*			Differential Evolution Methods, Mishra, SK (2007)
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:
#ifndef FKEANE_BUMP_H
#define FKEANE_BUMP_H

#include "BenchmarkFunction.h"


class FKeane_Bump : public BenchmarkFunction
{
    public:
		FKeane_Bump(ParamMap &v);
        FKeane_Bump(const int rId, const int rDimNumber, string& rName);
        virtual ~FKeane_Bump();
		bool isValid(const VirtualEncoding &ss);
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};

#endif // FKEANE_BUMP_H