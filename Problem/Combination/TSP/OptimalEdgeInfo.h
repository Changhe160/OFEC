/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com  Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Oct 2014
// Last modified:

#ifndef OptimalEdgeInfo_H
#define OptimalEdgeInfo_H

#include "../../../Global/global.h"
#include "../../../Global/solution.h"
#include "TravellingSalesman.h"
/* the ratio of the global optimal edges in the popualtion and the best individual during the search process and at the last iteration, and
the value of the last evaluation for TSP*/

class OptimalEdgeInfo
{
public:
	vector<vector<pair<int,double> > > mvv_edgeInfoPop;
	vector<vector<double> > mvv_edgeInfoBest;
	vector<pair<double,double> > mv_lastEdgeInfo;
	vector<int> mv_lastEvals;
	stringstream m_fileName;

	vector<vector<pair<int,int> > > mvv_diffEdges;
	vector<vector<double> > mvv_impRatio;
	void recordiffAndImp(int id,int eval,int diffs,double imp);

	OptimalEdgeInfo(ParamMap &v);
	~OptimalEdgeInfo(){}
	static void initialize(ParamMap &v);
	void setFileName(ParamMap &v);

	template<typename TypeIndiv>
	void recordEdgeInfo(Global *glob, const Solution<CodeVInt> & bestSofar,const vector<unique_ptr<TypeIndiv> > &pop, int &num,int popsize,int saveFre,bool flag=true)
	{
		if(glob->mp_problem->getEvaluations()/saveFre<num&&flag) return;

		vector<vector<double> > edgeInfo;
		int i,j;
		int numDim=glob->mp_problem->getNumDim();
		edgeInfo.resize(numDim);
		for(i=0;i<numDim;i++)
			edgeInfo[i].resize(numDim);

		double edgePop=0,edgeBest=0;
	
		for(i=0;i<numDim;i++)
			for(j=0;j<numDim;j++)
				edgeInfo[i][j]=0;
		
		for(i=0;i<numDim;i++)
		{
			edgeInfo[bestSofar.data().m_x[i]][bestSofar.data().m_x[(i+1)%numDim]]=1;
			edgeInfo[bestSofar.data().m_x[(i+1)%numDim]][bestSofar.data().m_x[i]]=1;
		}
		Solution<CodeVInt> & gopt=dynamic_cast<TravellingSalesman*>(glob->mp_problem.get())->getGOpt()[0];
		for(i=0;i<numDim;i++)
			if(edgeInfo[gopt.data().m_x[i]][gopt.data().m_x[(i+1)%numDim]])		edgeBest+=1;
		edgePop=bestSofar.getObjDistance_(gopt.data().m_obj)/gopt.data().m_obj[0];
		edgeBest/=numDim;

		edgePop=0;
		for(i=0;i<popsize;i++)
		{
			for(j=0;j<numDim;j++)
			{
				edgeInfo[int(pop[i]->data()[j])][int(pop[i]->data()[(j+1)%numDim])]=1;
				edgeInfo[int(pop[i]->data()[(j+1)%numDim])][int(pop[i]->data()[j])]=1;
			}
		}
		for(i=0;i<numDim;i++)
			if (edgeInfo[int(gopt.data().m_x[i])][int(gopt.data().m_x[(i + 1) % numDim])])
				edgePop+=1;
		edgePop/=numDim;


		mvv_edgeInfoPop[glob->m_runId].push_back(move(make_pair(glob->mp_problem->getEvaluations(),edgePop)));
		mvv_edgeInfoBest[glob->m_runId].push_back(edgeBest);
		++num;
	}

	template<typename TypeIndiv>
	void recordEdgeInfo(Global *glob, const vector<TypeIndiv> &pop, const vector<int>& bestIdx,int &num,int popsize,int saveFre,bool flag=true)
	{
		if(glob->mp_problem->getEvaluations()/saveFre<num&&flag) return;

		vector<vector<double> > edgeInfo;
		int i,j;
		int numDim=glob->mp_problem->getNumDim();
		edgeInfo.resize(numDim);
		for(i=0;i<numDim;i++)
			edgeInfo[i].resize(numDim);
		for(i=0;i<numDim;i++)
			for(j=0;j<numDim;j++)
				edgeInfo[i][j]=0;
		for(i=0;i<popsize;i++)
		{
			for(j=0;j<numDim;j++)
			{
				edgeInfo[int(pop[i].data()[j])][int(pop[i].data()[(j+1)%numDim])]=1;
				edgeInfo[int(pop[i].data()[(j+1)%numDim])][int(pop[i].data()[j])]=1;
			}
		}
		double edgePop=0,edgeBest=0;
		Solution<CodeVInt> & gopt=dynamic_cast<TravellingSalesman*>(glob->mp_problem.get())->getGOpt()[0];

		for(i=0;i<numDim;i++)
			if(edgeInfo[gopt.data().m_x[i]][gopt.data().m_x[(i+1)%numDim]])
				edgePop+=1;

		for(i=0;i<numDim;i++)
			for(j=0;j<numDim;j++)
				edgeInfo[i][j]=0;
		for(i=0;i<numDim;i++)
		{
			edgeInfo[int(pop[bestIdx[0]].data()[i])][int(pop[bestIdx[0]].data()[(i+1)%numDim])]=1;
			edgeInfo[int(pop[bestIdx[0]].data()[(i+1)%numDim])][int(pop[bestIdx[0]].data()[i])]=1;
		}
		for(i=0;i<numDim;i++)
			if(edgeInfo[gopt.data().m_x[i]][gopt.data().m_x[(i+1)%numDim]])
				edgeBest+=1;
		edgePop/=numDim;
		edgeBest/=numDim;
		mvv_edgeInfoPop[glob->m_runId].push_back(move(make_pair(glob->mp_problem->getEvaluations(),edgePop)));
		mvv_edgeInfoBest[glob->m_runId].push_back(edgeBest);
		++num;
	}

	void recordLastInfo(int ID, int Evals);
	void output();
	static OptimalEdgeInfo* getOptimalEdgeInfo();
	static void deleteOptimalEdgeInfo();
protected:
	static unique_ptr<OptimalEdgeInfo> msp_perf;
};

#endif