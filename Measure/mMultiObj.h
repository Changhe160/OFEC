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
// Created: 31 December 2014
// Last modified:

#ifndef MULTIOBJSTA_H
#define MULTIOBJSTA_H

#include "../Global/global.h"
#include "../Utility/myVector.h"

class mMultiObj
{
public:
	static mMultiObj* getMultiObj();
	static void deleteMultiObj();
	static void initialize(Global *,ParamMap &v);
	void record(Global * glob,int index,vector<double> &obj,vector<double> &point);
	template<typename TypeIndi>
	void recordDistance(Global *glob,int ID,const vector<unique_ptr<TypeIndi> > &pop)
	{
		int evals=glob->mp_problem->getEvaluations();
		if(evals%(int)Global::g_arg[param_sampleFre]==0)
		{
			double distance = 0;
			
			const vector<vector<double>*> &obj = glob->mp_problem->getObjSet();
			MyVector v2(glob->mp_problem->getNumObj());
			for(decltype(obj.size()) i=0; i<obj.size(); i++)
			{
				double min_d = 1.0e+10;		
				for(int j=0; j<pop.size(); j++)
				{
					v2=pop[j]->data().m_obj;
					double d = v2.getDis(*obj[i]);
					if(d<min_d)  min_d = d;
				}
				distance+= min_d;
			}
			distance/=obj.size();

			mvv_distance[ID].push_back(evals);
			mvv_distance[ID].push_back(distance);
		}
	}

	void recordEvals(Global *,int ID);
	void setFileName(ParamMap &v);
	void outputResult();
	void reInitialize(Global * glob,int pops);
	vector<vector<vector<double> > > mvvv_obj;
	vector<vector<vector<double> > > mvvv_point;
	vector<vector<double> > mvv_distance;
	double getCurDis2PF(int ridx) {
		return mvv_distance[ridx].back();
	}
	void recordSuccess(bool suc) {
		m_suc+=suc;
	}
protected:
	static unique_ptr<mMultiObj>  msp_perf;
	stringstream  m_fileName;
	mMultiObj(Global *,ParamMap &v);
	int m_suc=0;	// the number of runs that finds the POS 
	vector<int> m_psr;
	vector<vector<bool>> m_fpsr;
};

#endif //MULTIOBJSTA_H