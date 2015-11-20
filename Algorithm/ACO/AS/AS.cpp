/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Oct 2014
// Last modified:

#include<float.h>
#include "AS.h"
#include "../../../Global/global.h"
#include "../../../Problem/Combination/TSP/OptimalEdgeInfo.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
extern bool g_algTermination;
#endif

AS::AS(double alpha,double beta,double Q,int Popsize,int NC,int numDim,double coeff):m_alpha(alpha),m_beta(beta),m_Q(Q),\
	m_NC(NC),m_coeff(coeff),m_saveFre(1500),m_num(0),Population<CodeVInt,Ant>(Popsize),m_preAveVal(DBL_MAX),m_stopCount(0)
{
	mvv_phero.resize(numDim);
	for(int i=0;i<numDim;i++)
		mvv_phero[i].resize(numDim);
	m_stopCriterion=MIN_COVER;
}

AS::AS(ParamMap &v):m_alpha(1),m_beta(5),m_Q(100),m_NC(2500),m_coeff(0.5),m_saveFre(1500),m_num(0),Population<CodeVInt,Ant>(int(v[param_popSize])),\
	m_preAveVal(DBL_MAX),m_stopCount(0)
{
	mvv_phero.resize(int(v[param_numDim]));
	for(int i=0;i<int(v[param_numDim]);i++)
		mvv_phero[i].resize(int(v[param_numDim]));
	m_stopCriterion=MIN_COVER;
}

AS::~AS()
{
	mvv_phero.clear();
}

void AS::initializeSystem()
{
	int dim=m_pop[0]->getNumDim();
	for(int i=0;i<dim;i++)
		for(int j=0;j<dim;j++)
			mvv_phero[i][j]=1./dim;  
	for(int i=0;i<m_popsize;i++)
		m_pop[i]->initialization();
}

void AS::updatePheromeno()
{
	ReturnFlag rf;
	int i,j,dim;
	dim=m_pop[0]->getNumDim();
	for(i=0;i<m_popsize;i++)
	{
		rf=m_pop[i]->evaluate();
		if(rf==Return_Terminate) break;
	}
	#ifdef OFEC_CONSOLE
	vector<int> bestIdx=findBest();
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(),Solution<CodeVInt>::getBestSolutionSoFar(),m_pop,bestIdx,m_num,m_popsize,m_saveFre);
	#endif
	if(rf==Return_Terminate) return;
	vector<vector<double> > phero;
	phero.resize(dim);
	for(i=0;i<dim;i++)
		phero[i].resize(dim);
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			phero[i][j]=0;
	for(j=0;j<m_popsize;j++)
	{
		for(i=0;i<dim;i++)
		{
			phero[m_pop[j]->data()[i]][m_pop[j]->data()[(i+1)%dim]]+=m_Q/m_pop[j]->obj(0);
			phero[m_pop[j]->data()[(i+1)%dim]][m_pop[j]->data()[i]]=phero[m_pop[j]->data()[i]][m_pop[j]->data()[(i+1)%dim]]; //symmetric
		}
	}
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			mvv_phero[i][j]=m_coeff*mvv_phero[i][j]+phero[i][j];
}

void AS::resetAntsInfo()
{
	for(int i=0;i<m_popsize;i++)
		m_pop[i]->resetData();
}


bool AS::ifTerminating()
{
	#ifdef OFEC_DEMON
	if(g_algTermination) return true;
	#endif
	if(m_stopCriterion==MIN_COVER){
		if(m_stopCount>=200) return true;
	}else {
		if(m_evoNum>=m_NC) return true;
	}
	if(CAST_TSP->getGOpt().flagLoc())
	if(Solution<CodeVInt>::getBestSolutionSoFar().getObjDistance_(CAST_TSP->getGOpt()[0].data().m_obj)==0){
		if(m_evoNum==0) return false;
		return true;
	}
	
	return false;
}


ReturnFlag AS::run_()
{
	int i,j,dim;
	initializeSystem();
	dim=m_pop[0]->getNumDim();
	m_evoNum=0;
	while(!ifTerminating())
	{
		getStopCount();
		#ifdef OFEC_DEMON
			for(i=0;i<this->getPopSize();i++)
				updateBestArchive(this->m_pop[i]->self());
			vector<Algorithm*> vp;	
			vp.push_back(this);	
			msp_buffer->updateBuffer_(&vp);
		#endif
		for(i=1;i<dim;i++)
			for(j=0;j<m_popsize;j++)
				m_pop[j]->selectNextCity_Pro(mvv_phero,m_beta,m_alpha);
		updatePheromeno();
		resetAntsInfo();
		++m_evoNum;
		
		//cout<<" "<<Global::msp_global->mp_problem->getBestSolutionSoFar().getObjDistance(CAST_TSP->getGOpt()[0].data().m_obj)<<" "<<m_stopCount<<endl;
	}
	#ifdef OFEC_CONSOLE
	vector<int> bestIdx=findBest();
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(),Solution<CodeVInt>::getBestSolutionSoFar(),m_pop,bestIdx,m_num,m_popsize,m_saveFre,false);

	OptimalEdgeInfo::getOptimalEdgeInfo()->recordLastInfo(Global::msp_global->m_runId, Global::msp_global->mp_problem->getEvaluations());
	#endif
	return Return_Terminate;
}

void AS::getStopCount(){
	double currentAveValue=0;
	for(int i=0;i<m_popsize;i++)
		currentAveValue+=m_pop[i]->obj(0);
	currentAveValue/=m_popsize;
	if(m_preAveVal-currentAveValue>0.01)
	{	
		m_preAveVal=currentAveValue;	
		m_stopCount=0;
	}
	else m_stopCount++;
	
}