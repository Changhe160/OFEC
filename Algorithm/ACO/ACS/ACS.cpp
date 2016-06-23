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

Dorigo, M. (1996). "Ant Colony System: A Cooperative Learning Approach to the Traveling Salesman Problem"
IEEE Transactions on Evolutionary Computation.

*************************************************************************/
// Created: 1 Nov 2015
// Last modified:

#include<float.h>
#include "ACS.h"
#include "../../../Global/global.h"
#include "../../../Problem/Combination/TSP/OptimalEdgeInfo.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
extern bool g_algTermination;
#endif

ACS::ACS(double beta, double Q, int Popsize, int NC, int numDim, double coeff) :m_beta(beta), m_Q(Q), m_isHaveGlobalBest(false), \
m_NC(NC), m_coeff(coeff), m_saveFre(1500), m_num(0), Population<CodeVInt, Ant>(Popsize), m_globalBest(), m_impRadio(0)
{
	mvv_phero.resize(numDim);
	for (int i = 0; i<numDim; i++)
		mvv_phero[i].resize(numDim);
	m_stopCriterion = MIN_COVER;
	m_t = getLenOfNN();
	m_t = 1.0 / (numDim * m_t);

	if (m_stopCriterion == MIN_COVER){
		m_term.reset(new TermMean(Global::g_arg));
	}
	else if (m_stopCriterion == MAX_GEN){
		m_term.reset(new TermMaxGen(Global::g_arg));
	}
}

ACS::ACS(ParamMap &v) : m_beta(2), m_Q(0.9), m_NC(2500), m_coeff(0.1), m_saveFre(1500), m_num(0), Population<CodeVInt, Ant>(int(v[param_popSize])), \
m_globalBest(), m_isHaveGlobalBest(false), m_impRadio(0)
{
	mvv_phero.resize(int(v[param_numDim]));
	for (int i = 0; i<int(v[param_numDim]); i++)
		mvv_phero[i].resize(int(v[param_numDim]));
	m_stopCriterion = MIN_COVER;
	m_t = getLenOfNN();
	m_t = 1.0 / (int(v[param_numDim]) * m_t);

	if (m_stopCriterion == MIN_COVER){
		m_term.reset(new TermMean(v));
	}
	else if (m_stopCriterion == MAX_GEN){
		m_term.reset(new TermMaxGen(v));
	}
}

ACS::~ACS()
{
	mvv_phero.clear();
}

void ACS::initializeSystem()
{
	int dim = m_pop[0]->getNumDim();
	for (int i = 0; i<dim; i++)
		for (int j = 0; j<dim; j++)
			mvv_phero[i][j] = 1. / dim;
	for (int i = 0; i<m_popsize; i++)
		m_pop[i]->initialize(-1); 
}

void ACS::local_updatePheromeno(int ant_loc, bool isLastEdge)
{
	pair<int, int> edge;
	if (!isLastEdge)
		edge = m_pop[ant_loc]->getCurrentEdge();
	else 
		edge = m_pop[ant_loc]->getLastEdge();
	if (edge.first == edge.second && edge.first == -1)
		throw myException("edge error in local_updatePheromeno function @ACS.cpp");
	mvv_phero[edge.first][edge.second] = (1 - m_coeff)*mvv_phero[edge.first][edge.second] + m_coeff*m_t;
	mvv_phero[edge.second][edge.first] = (1 - m_coeff)*mvv_phero[edge.second][edge.first] + m_coeff*m_t;  //symmetric
}

void ACS::global_updatePheromeno() 
{
	ReturnFlag rf;
	int i, j, dim;
	dim = m_pop[0]->getNumDim();
	m_impRadio = 0;
	for (i = 0; i<m_popsize; i++)
	{
		double temp = m_pop[i]->data().m_obj[0];
		rf = m_pop[i]->evaluate();
		if (temp > m_pop[i]->data().m_obj[0])
			m_impRadio++;
		if (rf == Return_Terminate) break;
	}
	vector<int> bestIdx = findBest();
	if (!m_isHaveGlobalBest || m_globalBest < m_pop[bestIdx[0]]->representative()) //update globally best tour 
	{
		m_globalBest = m_pop[bestIdx[0]]->representative(); //deep copy
		m_isHaveGlobalBest = true;
	}
#ifdef OFEC_CONSOLE
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(), Solution<CodeVInt>::getBestSolutionSoFar(), m_pop, m_num, m_popsize, m_saveFre);
#endif
	if (rf == Return_Terminate) return;
	double len = 1./m_globalBest.data().m_obj[0];
	vector<vector<int>> global_edges(dim);
	for (int i = 0; i < dim; i++)
	{
		global_edges[i].resize(dim);
	}
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			global_edges[i][j] = 0;
	vector<int> edges(dim);
	for (int i = 0; i < dim; i++)
	{
		edges[i] = m_globalBest.data().m_x[i];
	}
	for (int i = 0; i < dim; i++)
	{
		global_edges[edges[i]][edges[(i+1)%dim]] = 1;
	}
	for (int i = 0; i < dim; i++)
		for (int j = i + 1; j < dim; j++) //symmetric
		{
			if (global_edges[i][j] || global_edges[j][i])
			{
				mvv_phero[i][j] = (1 - m_coeff)*mvv_phero[i][j] + m_coeff*len;
				mvv_phero[j][i] = mvv_phero[i][j];
			}
			else
			{
				mvv_phero[i][j] = (1 - m_coeff)*mvv_phero[i][j];
				mvv_phero[j][i] = mvv_phero[i][j];
			}
		}
}

void ACS::resetAntsInfo()
{
	for (int i = 0; i<m_popsize; i++)
		m_pop[i]->resetData();
}


bool ACS::ifTerminating()
{
	if (m_stopCriterion == MIN_COVER){
		if (dynamic_cast<TermMean*>(m_term.get())->ifTerminating()) {
			if (m_iter == 0) return false;
			return true;
		}
	}
	else if (m_stopCriterion == MAX_GEN){
		if (dynamic_cast<TermMaxGen*>(m_term.get())->ifTerminating(m_iter))
			return true;
	}
	return false;
}


ReturnFlag ACS::run_()
{
	int i, j, dim;
	initializeSystem();
	dim = m_pop[0]->getNumDim();
	m_iter = 0;

	if (m_stopCriterion == MIN_COVER){
		dynamic_cast<TermMean*>(m_term.get())->initialize(DBL_MAX);
	}

	while (!ifTerminating())
	{
#ifdef OFEC_DEMON
		for (i = 0; i<this->getPopSize(); i++)
			updateBestArchive(this->m_pop[i]->self());
		vector<Algorithm*> vp;
		vp.push_back(this);
		msp_buffer->updateBuffer_(&vp);
#endif
		for (i = 0; i < m_popsize; i++)
		{
			for (j = 1; j < dim; j++)
			{
				double q = Global::msp_global->mp_uniformAlg->Next();
				if (q <= m_Q)
					m_pop[i]->selectNextCity_Greedy(mvv_phero, m_beta);
				else
					m_pop[i]->selectNextCity_Pro(mvv_phero, m_beta);
				local_updatePheromeno(i);
			}
			local_updatePheromeno(i, true);
		}
		global_updatePheromeno();
		resetAntsInfo();
		++m_iter;
		if (m_stopCriterion == MIN_COVER) {
			dynamic_cast<TermMean*>(m_term.get())->countSucIter(mean());
		}
		//cout<<" "<<Global::msp_global->mp_problem->getBestSolutionSoFar().getObjDistance(CAST_TSP->getGOpt()[0].data().m_obj)<<" "<<m_stopCount<<endl;
#ifdef OFEC_CONSOLE
		double tempdif = 0;
		for (int i = 0; i < m_popsize; i++)
			tempdif += m_pop[i]->self().getDistance(Solution<CodeVInt>::getBestSolutionSoFar());
		tempdif /= m_popsize;
		double impr = static_cast<double>(m_impRadio) / m_popsize;
		OptimalEdgeInfo::getOptimalEdgeInfo()->recordiffAndImp(Global::msp_global->m_runId, Global::msp_global->mp_problem->getEvaluations(), fabs(tempdif), impr);
#endif
	}
#ifdef OFEC_CONSOLE
	vector<int> bestIdx = findBest();
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(), Solution<CodeVInt>::getBestSolutionSoFar(), m_pop, m_num, m_popsize, m_saveFre, false);

	OptimalEdgeInfo::getOptimalEdgeInfo()->recordLastInfo(Global::msp_global->m_runId, Global::msp_global->mp_problem->getEvaluations());
#endif
	return Return_Terminate;
}


double ACS::getLenOfNN()
{
	vector<int> candidate(m_numDim), result(m_numDim);
	TravellingSalesman *_ptr = dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get());
	const vector<vector<double>> cost = _ptr->getCost();
	int n = 0;
	for (int i = 0; i < candidate.size(); i++){
		candidate[i] = i;
	}
	result[n++] = candidate[0];
	candidate[0] = candidate[m_numDim - 1];
	while (n < m_numDim){
		int loc = 0;
		double min = cost[result[n - 1]][candidate[loc]];
		for (int m = 1; m < m_numDim - n; m++){
			if (cost[result[n - 1]][candidate[m]] < min){
				min = cost[result[n - 1]][candidate[m]];
				loc = m;
			}
		}
		result[n++] = candidate[loc];
		candidate[loc] = candidate[m_numDim - n];
	}
	double val = 0;
	for (int i = 0; i < m_numDim; i++){
		val += cost[result[i]][result[(i + 1) % m_numDim]];
	}
	return val;
}