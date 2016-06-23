#include "MMAS.h"
#include "../../../Global/global.h"
#include "../../../Problem/Combination/TSP/OptimalEdgeInfo.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
extern bool g_algTermination;
#endif

static mutex m_nearMutex;
static vector<set<int>> mv_candidate; //select a next city from the candidate list

MMAS::MMAS(double alpha, double beta, double Q, int Popsize, int NC, int numDim, double coeff) : AS(alpha, beta, Q, Popsize, NC, numDim, coeff),\
m_isHaveGlobalBest(false), m_isHaveRestartBest(false), m_globalBest(), m_restartBest(), m_impRadio(0)
{
	setDefaultParameters();
	if (m_globalBest.getNumDim() <= 20)
		m_length = m_globalBest.getNumDim() - 1;
	findnearghbor();
}

MMAS::MMAS(ParamMap &v) : AS(v), m_isHaveGlobalBest(false), m_isHaveRestartBest(false), m_globalBest(), m_restartBest(), m_impRadio(0)
{
	setDefaultParameters();
	if (m_globalBest.getNumDim() <= 20)
		m_length = m_globalBest.getNumDim() - 1;
	findnearghbor();
}

void MMAS::setDefaultParameters()
{
	m_alpha = 1.0;
	m_beta = 2.0;
	m_length = 20;
	m_rho = 0.02;
	m_branchFrc = 1.00001;
	m_lambda = 0.05;
	m_uGB = LONG_MAX;
	m_iter = 1;
	m_restartFoundBest = 0;

	TravellingSalesman * ptr = dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get());
	double temp = getLenOfNN();
	m_pheroMax = 1. / (m_rho * temp);
	m_pheroMin = m_pheroMax / (2. * m_globalBest.getNumDim());
}

double MMAS::nodeBranching()
{
	int n = m_globalBest.getNumDim();
	double minTemp, maxTemp, cutoff, avg = 0.;
	vector<double> num_branches(n);
	set<int>::iterator iter;
	for (int i = 0; i < n; i++)
	{
		iter = mv_candidate[i].begin();
		minTemp = mvv_phero[i][*iter];
		maxTemp = mvv_phero[i][*iter];
		++iter;
		for (; iter != mv_candidate[i].end(); iter++)
		{
			if (mvv_phero[i][*iter] > maxTemp)
				maxTemp = mvv_phero[i][*iter];
			if (mvv_phero[i][*iter] < minTemp)
				minTemp = mvv_phero[i][*iter];
		}
		cutoff = minTemp + m_lambda * (maxTemp - minTemp);
		num_branches[i] = 0.;
		for (set<int>::iterator it = mv_candidate[i].begin(); it != mv_candidate[i].end(); it++)
		{
			if (mvv_phero[i][*it] > cutoff)
				num_branches[i] += 1.;
		}
		avg += num_branches[i];
	}
	return avg / (n * 2.);
}

MMAS::~MMAS()
{

}

void MMAS::initializeSystem(double t0, bool isDefault)
{
	if (isDefault)
		AS::initializeSystem();
	else
		AS::initializeSystem(t0, false);
}

void MMAS::updatePheroMinAndMax()
{
	int n = m_globalBest.getNumDim();
	double p_x = exp(log(0.05) / n);
	m_pheroMin = 1. * (1. - p_x) / (p_x * (double)((m_length + 1) / 2));
	m_pheroMax = 1. / ((m_rho)* m_globalBest.data().m_obj[0]);
	m_pheroMin = m_pheroMax * m_pheroMin;
}

void MMAS::updatePheromeno()
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
		m_restartBest = m_pop[bestIdx[0]]->representative(); //deep copy
		m_isHaveGlobalBest = true;
		m_isHaveRestartBest = true;

		m_restartFoundBest = m_iter;
		m_branchFactor = nodeBranching();
		updatePheroMinAndMax();
	}
	if (!m_isHaveRestartBest || m_restartBest < m_pop[bestIdx[0]]->representative())
	{
		m_restartBest = m_pop[bestIdx[0]]->representative(); //deep copy
		m_restartFoundBest = m_iter;
		m_isHaveRestartBest = true;
	}
#ifdef OFEC_CONSOLE
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(), Solution<CodeVInt>::getBestSolutionSoFar(), m_pop, m_num, m_popsize, m_saveFre);
#endif
	if (rf == Return_Terminate) return;
	for (i = 0; i < dim; i++)
		for (j = i+1; j < dim; j++)
		{
			mvv_phero[i][j] = (1 - m_rho) * mvv_phero[i][j];
			mvv_phero[j][i] = mvv_phero[i][j];
		}
	if (m_iter % m_uGB)
	{
		double len = 1.0 / m_pop[bestIdx[0]]->data().m_obj[0];
		for (i = 0; i < dim; i++)
		{
			int cur = m_pop[bestIdx[0]]->data().m_x[i];
			int next = m_pop[bestIdx[0]]->data().m_x[(i + 1) % dim];
			mvv_phero[cur][next] += len;
			mvv_phero[next][cur] = mvv_phero[cur][next];
		}
	}
	else
	{
		double len = 1.0 / m_restartBest.data().m_obj[0];
		for (i = 0; i < dim; i++)
		{
			int cur = m_restartBest.data().m_x[i];
			int next = m_restartBest.data().m_x[(i + 1) % dim];
			mvv_phero[cur][next] += len;
			mvv_phero[next][cur] = mvv_phero[cur][next];
		}
	}
	m_uGB = 25;
	checkPheromoneTrailLimits();
}

void MMAS::checkPheromoneTrailLimits()
{
	int n = m_globalBest.getNumDim();

	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < i; j++)
		{
			if (mvv_phero[i][j] < m_pheroMin) 
			{
				mvv_phero[i][j] = m_pheroMin;
				mvv_phero[j][i] = m_pheroMin;
			}
			else if (mvv_phero[i][j] > m_pheroMax) {
				mvv_phero[i][j] = m_pheroMax;
				mvv_phero[j][i] = m_pheroMax;
			}
		}
	}
}

ReturnFlag MMAS::run_()  
{
	int i, j, dim;
	initializeSystem(m_pheroMax,false);
	dim = m_pop[0]->getNumDim();

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
		for (i = 1; i < dim; i++)
			for (j = 0; j < m_popsize; j++)
				m_pop[j]->selectNextCity_Pro(mvv_phero, mv_candidate, m_beta, m_alpha);
		updatePheromeno();
		pheroSmoothing();
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
	OptimalEdgeInfo::getOptimalEdgeInfo()->recordEdgeInfo<Ant>(Global::msp_global.get(), Solution<CodeVInt>::getBestSolutionSoFar(), m_pop, m_num, m_popsize, m_saveFre, false);

	OptimalEdgeInfo::getOptimalEdgeInfo()->recordLastInfo(Global::msp_global->m_runId, Global::msp_global->mp_problem->getEvaluations());
#endif
	return Return_Terminate;
}

void MMAS::findnearghbor()
{
	m_nearMutex.lock();

	if (mv_candidate.empty())
	{
		vector<vector<int>> temp;
		mv_candidate.resize(m_globalBest.getNumDim());
		temp.resize(m_globalBest.getNumDim());
		for (int i = 0; i < temp.size(); i++)
			temp[i].resize(m_length);
		dynamic_cast<TravellingSalesman*>(Global::msp_global->mp_problem.get())->findNearbyCity(temp);
		for (int i = 0; i < temp.size(); i++)
		{
			for (int j = 0; j < m_length; j++)
			{
				mv_candidate[i].insert(temp[i][j]);
			}
		}
	}

	m_nearMutex.unlock();
}

void MMAS::pheroSmoothing() 
{
	if (!(m_iter % 100))
	{
		m_branchFactor = nodeBranching();
		if (m_branchFactor < m_branchFrc && (m_iter - m_restartFoundBest > 250))
		{
			m_isHaveRestartBest = false;
			initPhero();
		}
	}
}

void MMAS::initPhero()
{
	int n = m_globalBest.getNumDim();
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			mvv_phero[i][j] = m_pheroMax;
			mvv_phero[j][i] = mvv_phero[i][j];
		}
	}
}

double MMAS::getLenOfNN()
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