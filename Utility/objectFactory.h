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
// Created: 21 September 2014
// Last modified: 12 Dec. 2014

#ifndef OBJECT_FACTORY_H
#define OBJECT_FACTORY_H

#include "include.h"
#include "TypeVar/typeVar.h"

class Problem;
class Algorithm;

struct classFactory{
	typedef Problem * (*constructorProblem)(ParamMap&);
	typedef map<string, constructorProblem> ClassMapProblem;

	typedef Algorithm * (*constructorAlgorithm)(ParamMap&);
	typedef map<string, constructorAlgorithm> ClassMapAlgorithm;

	constructorProblem constructProblem(const string &s);
	constructorAlgorithm constructAlgorithm(const string &s);

	ClassMapProblem m_theMapProblem;
	ClassMapAlgorithm m_theMapAlgorithm;
};



#endif