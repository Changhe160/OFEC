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
#include "mSingleObjDyn.h"
#include "../Global/global.h"

mSingleObjDyn::mSingleObjDyn(ParamMap &v):mSingleObj(v)
{
	m_recordsPerChange=m_changeFre/Global::g_arg[param_sampleFre];
}

void mSingleObjDyn::initialize(ParamMap &v){

    if(mSingleObj::msp_perf)        return;

    mSingleObj::msp_perf=unique_ptr<mSingleObj>(new mSingleObjDyn(v));
}

void mSingleObjDyn::calculateKeyParam()
{
	calculateNumRecords();
	int maxRun=MAX_NUM_RUN;
	m_numChanges=Global::g_arg[param_numChange];
	for(int i=0;maxRun>i;i++) mpp_convergeTime[i].resize(m_numChanges);
}