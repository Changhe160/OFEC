#ifndef MATH_AUX__
#define MATH_AUX__

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 11 Jan 2015
// Last modified:

#include <cstdlib>
#include <vector>

namespace MathAux
{
// ASF(): achievement scalarization function
double ASF(const std::vector<double> &objs, const std::vector<double> &weight);

// GuassianElimination(): used to calculate the hyperplane
void GuassianElimination(std::vector<double> *px, std::vector< std::vector<double> > A, const std::vector<double> &b);

// PerpendicularDistance(): calculate the perpendicular distance from a point to a line
double PerpendicularDistance(const std::vector<double> &direction, const std::vector<double> &point);
}

#endif