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

/**************************************************************************
design of a gear train,  which
was introduced in
    E.~Sandgren, ``Nonlinear integer and discrete programming in mechanical
    design,'' in \emph{the ASME Design Technology Conf.,}, 1988, pp. 95--105.
, is to optimize the gear ratio for a compound gear train that contains three
gears. It is to be designed that the gear ratio is as close as possible to
1/6.931. For each gear, the number of teeth must be between
12 and 60.
******************************************************************************/

// Created: 21 July 2011
// Last modified:
#ifndef FGEAR_TRAIN_H
#define FGEAR_TRAIN_H

#include "../FunctionOpt/BenchmarkFunction.h"

class FGear_Train : public BenchmarkFunction
{
    public:
		FGear_Train(ParamMap &v);
        FGear_Train(const int rId,  const int rDim, string& rName);
        virtual ~FGear_Train();
    protected:
        void initialize();
        void evaluate__(double const *x,vector<double>& obj);
    private:
};
#endif // FGEAR_TRAIN_H
