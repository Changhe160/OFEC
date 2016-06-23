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

#include "OptimalEdgeInfo.h"
#include "../../../Global/global.h"
#include "../../../Measure/mSingleObj.h"
#include "../TSP/TravellingSalesman.h"
#include "../../../Algorithm/Population.h"

unique_ptr<OptimalEdgeInfo> OptimalEdgeInfo::msp_perf(nullptr);
void OptimalEdgeInfo::initialize(ParamMap &v)
{
	if(OptimalEdgeInfo::msp_perf)        return;
    OptimalEdgeInfo::msp_perf=unique_ptr<OptimalEdgeInfo>(new OptimalEdgeInfo(v));
}

OptimalEdgeInfo::OptimalEdgeInfo(ParamMap &v)
{
	mvv_edgeInfoBest.resize(int(MAX_NUM_RUN));
	mvv_edgeInfoPop.resize(int(MAX_NUM_RUN));
	mv_lastEdgeInfo.resize(int(MAX_NUM_RUN));
	mv_lastEvals.resize(int(MAX_NUM_RUN));
	mvv_diffEdges.resize(int(MAX_NUM_RUN));
	mvv_impRatio.resize(int(MAX_NUM_RUN));
	setFileName(v);
}

void OptimalEdgeInfo::output()
{
	size_t dim=mvv_edgeInfoBest[0].size();
	for(int i=1;i<int(MAX_NUM_RUN);i++)
	{
		if(mvv_edgeInfoBest[i].size()>dim)
			dim=mvv_edgeInfoBest[i].size();
	}
	for(int i=0;i<int(MAX_NUM_RUN);i++)
	{
		int temp=mvv_edgeInfoBest[i].size();
		double edge=mvv_edgeInfoBest[i].back();
		pair<int,double> tempair=mvv_edgeInfoPop[i].back();
		for(int j=temp;j<dim;j++)
		{
			mvv_edgeInfoBest[i].push_back(edge);
			mvv_edgeInfoPop[i].push_back(tempair);
		}
	}
	vector<pair<double,double> > edgeInfoPop(dim);
	vector<double> edgeInfoBest(dim);
	pair<long double,double> tempPop;
	double tempBest;
	for(size_t i=0;i<dim;i++)
	{ 
		tempPop.first=0; tempPop.second=0;
		tempBest=0;
		for(int j=0;j<int(MAX_NUM_RUN);j++)
		{
			tempPop.first+=mvv_edgeInfoPop[j][i].first;
			tempPop.second+=mvv_edgeInfoPop[j][i].second;
			tempBest+=mvv_edgeInfoBest[j][i];
		}
		tempPop.first/=int(MAX_NUM_RUN); tempPop.second/=int(MAX_NUM_RUN);
		edgeInfoPop[i]=tempPop;
		edgeInfoBest[i]=tempBest/int(MAX_NUM_RUN);
	}
	ostringstream os1,os2,os3,os4;
	ofstream out1,out2,out3,out4;
	os1<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Optimumedges.txt";
	os2<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"lastOptimaEdgePro.txt";
	os3<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"lastObjective.txt";
	os4<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Sta.txt";
	out1.open(os1.str().c_str());
	out2.open(os2.str().c_str());
	out3.open(os3.str().c_str());
	out4.open(os4.str().c_str(),ios::app);
	if(!out1)
		throw myException("error in @OptimalEdgeInfo::output()");
	double lastPop=0,lastBest=0;
	double lastError=0,lastEval=0;
	for(int i=0;i<int(MAX_NUM_RUN);i++)
	{
		lastPop+=mv_lastEdgeInfo[i].first;
		lastBest+=mv_lastEdgeInfo[i].second;
		lastError+=(mSingleObj::getSingleObj()->getBestSoFar(i)-mSingleObj::getSingleObj()->getGOpt()[0][0]);
		lastEval+=mv_lastEvals[i];
	}
	lastPop/=int(MAX_NUM_RUN);
	lastBest/=int(MAX_NUM_RUN);
	lastError/=int(MAX_NUM_RUN);
	lastEval/=int(MAX_NUM_RUN);
	for(size_t i=0;i<dim;i++)
		out1<<edgeInfoPop[i].second<<" "<<edgeInfoBest[i]<<" "<<edgeInfoPop[i].first<<endl;
	for(int i=0;i<int(MAX_NUM_RUN);i++)
	{
		out2<<mv_lastEdgeInfo[i].first<<" "<<mv_lastEdgeInfo[i].second<<endl;
		out3<<mSingleObj::getSingleObj()->getBestSoFar(i)-mSingleObj::getSingleObj()->getGOpt()[0][0]<<" "<<mv_lastEvals[i]<<endl;
	}
	out2<<"ave: "<<lastPop<<" "<<lastBest;
	out3<<"ave: "<<lastError<<" "<<lastEval;
	out4<<"Mean of edgeInfoPop: "<<lastPop<<endl;
	out4<<"Mean of edgeInfoBest: "<<lastBest<<endl;
	out4<<"Mean of lastEvals: "<<lastEval<<endl;

	out1.close(); out1.clear();
	out2.close(); out2.clear();
	out3.close(); out3.clear();
	out4.close(); out4.clear();

	int max=mvv_diffEdges[0].size();
	for(int i=1;i<int(MAX_NUM_RUN);i++)
		if(mvv_diffEdges[i].size()>max)
			max=mvv_diffEdges[i].size();
	for(int i=0;i<int(MAX_NUM_RUN);i++)
	{
		if(mvv_diffEdges[i].size()<max)
		{
			pair<int,int> last=mvv_diffEdges[i].back();
			for(int j=mvv_diffEdges[i].size();j<max;j++)
				mvv_diffEdges[i].push_back(last);
		}
		if(mvv_impRatio[i].size()<max)
		{
			double last=mvv_impRatio[i].back();
			for(int j=mvv_impRatio[i].size();j<max;j++)
				mvv_impRatio[i].push_back(last);
		}
	}
	vector<double> impRatio_ave(max,0);
	vector<pair<double,double> > diffEdges(max,make_pair(0,0));
	for(int j=0;j<max;j++)
	{
		for(int i=0;i<int(MAX_NUM_RUN);i++)
		{
			impRatio_ave[j]+=mvv_impRatio[i][j];
			diffEdges[j].first+=mvv_diffEdges[i][j].first;
			diffEdges[j].second+=mvv_diffEdges[i][j].second;
		}
		impRatio_ave[j]/=int(MAX_NUM_RUN);
		diffEdges[j].first/=int(MAX_NUM_RUN);
		diffEdges[j].second/=int(MAX_NUM_RUN);
	}
	os1.str("");
	os1<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"ImpRatioAndDiffs.txt";
	out1.open(os1.str().c_str());
	for(int i=0;i<max;i++)
		out1<<impRatio_ave[i]<<" "<<diffEdges[i].first<<" "<<diffEdges[i].second<<endl;
	out1.close();
}

void OptimalEdgeInfo::setFileName(ParamMap &v){
	m_fileName.str("");
	
	for(auto &i:v){
		for(auto &j:Global::msm_param){
			if(i.first==param_gOptFlag||i.first==param_algId||i.first==param_proId||i.first==param_flagNoise||\
				i.first==param_flagNumPeakChange||i.first==param_flagTimeLinkage||i.first==param_numRun||\
				i.first==param_numTask||i.first==param_minNumPopSize||i.first==param_hibernatingRadius||\
				i.first==param_solutionValidationMode||i.first==param_evalCountFlag||\
				i.first==param_workingDir||i.first==param_sampleFre||i.first==param_maxEvals) continue;
			if(i.first==j.second){			
				m_fileName<<j.first.substr(6)<<i.second<<"_";			
				break;
			}
		}		
	}
}


void OptimalEdgeInfo::recordLastInfo(int ID, int Evals)
{
	mvv_edgeInfoPop[ID].push_back(make_pair(Evals,mvv_edgeInfoPop[ID].back().second));
	mvv_edgeInfoBest[ID].push_back(mvv_edgeInfoBest[ID].back());
	mv_lastEdgeInfo[ID]=make_pair(mvv_edgeInfoPop[ID].back().second,mvv_edgeInfoBest[ID].back());
	mv_lastEvals[ID]=Evals;
}

void OptimalEdgeInfo::recordiffAndImp(int id,int eval,int diffs,double imp)
{
	mvv_impRatio[id].push_back(imp);
	mvv_diffEdges[id].push_back(make_pair(diffs,eval));
}


OptimalEdgeInfo* OptimalEdgeInfo::getOptimalEdgeInfo(){
	 return OptimalEdgeInfo::msp_perf.get();
}

void OptimalEdgeInfo::deleteOptimalEdgeInfo(){
	return OptimalEdgeInfo::msp_perf.reset();
}

