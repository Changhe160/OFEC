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
// Created: 21 July 2011
// Last modified:

#ifndef SINGLEOBJDYN_H
#define SINGLEOBJDYN_H

#include "mSingleObj.h"

class mSingleObjDyn: public mSingleObj
{
public:
	static void initialize(ParamMap &v);
	mSingleObjDyn(ParamMap &v);
	~mSingleObjDyn(){}
protected:
	void calculateKeyParam();
};


#endif //SINGLEOBJDYN_H