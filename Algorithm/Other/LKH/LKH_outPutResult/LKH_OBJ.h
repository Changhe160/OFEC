/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 Jan 2015
// Last modified:

#ifndef MLKHOBJ_H
#define MLKHOBJ_H

#include "../LKH.h"

class mLKHObj 
{
public:
	static void initialize(ParamMap &v);
	static mLKHObj* getLKHOBJ();
	static void deleteLKHOBJ();
	void outputResult();
	stringstream m_fileName;
protected:
	void setFileName(ParamMap &v);
	mLKHObj(ParamMap &v);
	static unique_ptr<mLKHObj> msp_perf;
};

#endif //MLKHOBJ_H