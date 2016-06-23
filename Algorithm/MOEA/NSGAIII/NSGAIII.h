#ifndef NSGAIII__
#define NSGAIII__

#include <cstddef>
#include <string>
#include <fstream>
#include <vector>
#include "ReferencePoint.h"
#include "MathAux.h"
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
// ----------------------------------------------------------------------------------

//		NSGAIII
//
// Deb and Jain, "An Evolutionary Many-Objective Optimization Algorithm Using 
// Reference-point Based Non-dominated Sorting Approach, Part I: Solving Problems with 
// Box Constraints," IEEE Transactions on Evolutionary Computation, to appear.
//
// http://dx.doi.org/10.1109/TEVC.2013.2281535
// ----------------------------------------------------------------------------------
*************************************************************************/
// Created: 11 Jan 2015
// Last modified:

#include "../../Population.h"
#include "../../Individual.h"
#include "../../../Measure/mMultiObj.h"
#include "../../../Problem/FunctionOpt/MOP/DTLZ/DTLZ.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

template<typename TypeIndiv,typename TypePop>
class NSGAIII : public Algorithm
{
public:
	NSGAIII(ParamMap &v);
	~NSGAIII();
	ReturnFlag run_();
	virtual void evolve_mo()=0;

protected:
	void EnvironmentalSelection();
	vector<double> TranslateObjectives(const vector<vector<int> > &fronts);
	void FindExtremePoints(vector<size_t> *extreme_points, const vector<vector<int> > &fronts);
	void ConstructHyperplane(vector<double> *pintercepts, const vector<size_t> &extreme_points);
	void NormalizeObjectives(const vector<vector<int> > &fronts, const vector<double> &intercepts, const vector<double> &ideal_point);
	size_t FindNicheReferencePoint(const vector<ReferencePoint> &rps);
	int SelectClusterMember(const ReferencePoint &rp);
	void setdefaultparam();

	vector<size_t> mv_obj_division_p_; 
	vector<ReferencePoint> mv_rps;
	TypePop m_parent, m_offspring;
	vector<vector<double> > mvv_offConvObj;

};

template<typename TypeIndiv,typename TypePop>
NSGAIII<TypeIndiv,TypePop>::NSGAIII(ParamMap &v) : Algorithm(-1,string()), m_parent(),m_offspring(),mv_rps(0)
{
	int numObj=Global::msp_global->mp_problem->getNumObj();
	setdefaultparam();
	GenerateReferencePoints(&mv_rps, numObj, mv_obj_division_p_); 
	size_t PopSize = mv_rps.size();
	while (PopSize%4) PopSize += 1;

	if(m_parent.getPopSize()>0)
		m_parent.remove(m_parent.getPopSize());
	if(m_offspring.getPopSize()>0)
		m_offspring.remove(m_parent.getPopSize());

	TypePop temp1(PopSize,true);
	TypePop temp2(2*PopSize,false);
	m_parent.add(temp1);
	m_offspring.add(temp2);
	mvv_offConvObj.resize(2*PopSize);
	for(int i=0;i<2*PopSize;i++)
		mvv_offConvObj[i].resize(numObj);

	m_parent.clearBestArchive();
	m_offspring.clearBestArchive();

}

template<typename TypeIndiv,typename TypePop>
NSGAIII<TypeIndiv,TypePop>::~NSGAIII()
{
	mv_rps.clear();
	mvv_offConvObj.clear();
}

template<typename TypeIndiv,typename TypePop>
void NSGAIII<TypeIndiv,TypePop>::setdefaultparam()
{
	int numObj=Global::msp_global->mp_problem->getNumObj();
	if(gGetProblemName(Global::ms_curProId).find("FUN_MOP_DTLZ")!=string::npos)
	{
		if(numObj==3) 
			mv_obj_division_p_.resize(1,12);
		else if(numObj==5)
			mv_obj_division_p_.resize(1,6);
		else if(numObj==8)
		{
			mv_obj_division_p_.resize(2);
			mv_obj_division_p_[0]=3;
			mv_obj_division_p_[1]=2;
		}
		return ;
	}
	mv_obj_division_p_.resize(1,12);
}

template<typename TypeIndiv,typename TypePop>
ReturnFlag NSGAIII<TypeIndiv,TypePop>::run_()
{
#ifdef OFEC_CONSOLE
	if(mMultiObj::getMultiObj()&&Global::msp_global->mp_problem->isGlobalOptKnown())
		mMultiObj::getMultiObj()->recordDistance<TypeIndiv>(Global::msp_global.get(),Global::msp_global->m_runId,m_parent.getPop());
#endif

	// evolution
	while(!this->ifTerminating())
	{
		//cout << "Run " << Global::msp_global->m_runId << "  " << Global::msp_global->mp_problem->getEvaluations() << " " << mMultiObj::getMultiObj()->getCurDis2PF(Global::msp_global->m_runId) << endl;
		
		evolve_mo();
		EnvironmentalSelection();
		
#ifdef OFEC_DEMON
		vector<Algorithm*> vp;
		this->m_parent.rank();
		vp.push_back(&this->m_parent);
		msp_buffer->updateBuffer_(&vp);
#endif	

		#ifdef OFEC_CONSOLE
		if(mMultiObj::getMultiObj()&&Global::msp_global->mp_problem->isGlobalOptKnown())
			mMultiObj::getMultiObj()->recordDistance<TypeIndiv>(Global::msp_global.get(),Global::msp_global->m_runId,m_parent.getPop());
#endif

	}
	#ifdef OFEC_CONSOLE
	if(mMultiObj::getMultiObj()){
	mMultiObj::getMultiObj()->reInitialize(Global::msp_global.get(),m_parent.getPopSize());
	for(int i=0;i<m_parent.getPopSize();i++)
		mMultiObj::getMultiObj()->record(Global::msp_global.get(),i,m_parent[i]->data().m_obj,m_parent[i]->data().m_x);
	}
#endif

//	cout<<"Run "<<Global::msp_global->m_runId<<" is terminated"<<endl;
	return Return_Normal;
}

// ----------------------------------------------------------------------
// TranslateObjectives():
//
// 1. Find the ideal point
// 2. Translate the objective values
// 3. Return the ideal point
//
// Check steps 1-3 in Algorithm 2 in the original paper of NSGAIII.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
vector<double> NSGAIII<TypeIndiv,TypePop>::TranslateObjectives(const vector<vector<int> > &fronts)
{
	int numObj=Global::msp_global->mp_problem->getNumObj();
	vector<double> ideal_point(numObj);

	for (int f=0; f<numObj; f+=1)
	{
		double minf = numeric_limits<double>::max();
		for (size_t i=0; i<fronts[0].size(); i+=1) // min values must appear in the first front
		{
			minf = std::min(minf, m_offspring[ fronts[0][i] ]->data().m_obj[f]);
		}
		ideal_point[f] = minf;

		for (size_t t=0; t<fronts.size(); t+=1)
		{
			for (size_t i=0; i<fronts[t].size(); i+=1)
			{
				size_t ind = fronts[t][i];
				mvv_offConvObj[ind][f] = m_offspring[ind]->data().m_obj[f] - minf;
			}
		}
	}

	return ideal_point;

}// TranslateObjectives()

// ----------------------------------------------------------------------
// FindExtremePoints():
// 
// Find the extreme points along each objective axis.
// The extreme point has the minimal ASF value.
// Return the indices of extreme individuals in the population.
//
// Check step 4 in Algorithm 2 and eq. (4) in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
void NSGAIII<TypeIndiv,TypePop>::FindExtremePoints(vector<size_t> *extreme_points, const vector<vector<int> > &fronts)
{
	int numObj=Global::msp_global->mp_problem->getNumObj();
	vector<size_t> &exp = *extreme_points;
	exp.clear();

	for (size_t f=0; f<numObj; f+=1)
	{
		vector<double> w(numObj, 0.000001);
		w[f] = 1.0;

		double min_ASF = numeric_limits<double>::max();
		size_t min_indv = fronts[0].size();

		for (size_t i=0; i<fronts[0].size(); i+=1)  // only consider the individuals in the first front
		{
			double asf = MathAux::ASF(m_offspring[ fronts[0][i] ]->data().m_obj, w);
			if ( asf < min_ASF )
			{
				min_ASF = asf;
				min_indv = fronts[0][i];
			}
		}
		
		exp.push_back(min_indv);
	}

}// FindExtremePoints()

// ----------------------------------------------------------------------
// ConstructHyperplane():
//
// Given the extreme points, construct the hyperplane.
// Then, calculate the intercepts.
//
// Check step 6 in Algorithm 2 in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
void NSGAIII<TypeIndiv,TypePop>::ConstructHyperplane(vector<double> *pintercepts, const vector<size_t> &extreme_points)
{
	// Check whether there are duplicate extreme points.
	// This might happen but the original paper does not mention how to deal with it.
	int numObj=Global::msp_global->mp_problem->getNumObj();
	bool duplicate = false;
	for (size_t i=0; !duplicate && i<extreme_points.size(); i+=1)
	{
		for (size_t j=i+1; !duplicate && j<extreme_points.size(); j+=1)
		{
			duplicate = (extreme_points[i] == extreme_points[j]);
		}
	}

	vector<double> &intercepts = *pintercepts;
	intercepts.assign(numObj, 0);

	if (duplicate) // cannot construct the unique hyperplane (this is a casual method to deal with the condition)
	{
		for (size_t f=0; f<intercepts.size(); f+=1)
		{
			// extreme_points[f] stands for the individual with the largest value of objective f
			intercepts[f] = m_offspring[ extreme_points[f] ]->data().m_obj[f]; 
		}
	}
	else
	{
		// Find the equation of the hyperplane
		vector<double> b(numObj, 1.0);
		vector< vector<double> > A;
		for (size_t p=0; p<extreme_points.size(); p+=1)
		{
			A.push_back(m_offspring[ extreme_points[p] ]->data().m_obj);
		}
		vector<double> x;
		MathAux::GuassianElimination(&x, A, b);
	
		// Find intercepts
		for (size_t f=0; f<intercepts.size(); f+=1)
		{
			intercepts[f] = 1.0/x[f];
		}
	}
}

// ----------------------------------------------------------------------
// NormalizeObjectives():
// 
// Normalize objective values with respect to the intercepts and the ideal point.
// Check step  7 in Algorithm 2 and eq. (5) in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
void NSGAIII<TypeIndiv,TypePop>::NormalizeObjectives(const vector<vector<int> > &fronts, const vector<double> &intercepts, const vector<double> &ideal_point)
{	
	int numObj=Global::msp_global->mp_problem->getNumObj();
	for (size_t t=0; t<fronts.size(); t+=1)
	{
		for (size_t i=0; i<fronts[t].size(); i+=1)
		{
			size_t ind = fronts[t][i];
			for (size_t f=0; f<numObj; f+=1)
			{
				if ( fabs(intercepts[f] - ideal_point[f])>10e-10 ) // avoid the divide-by-zero error
					mvv_offConvObj[ ind ][f] = mvv_offConvObj[ ind ][f]/(intercepts[f] - ideal_point[f]);
				else
					mvv_offConvObj[ ind ][f] = mvv_offConvObj[ ind ][f]/10e-10;
			}
		}
	}

}// NormalizeObjectives()

// ----------------------------------------------------------------------
// FindNicheReferencePoint():
// 
// Find the reference point with the minimal cluster size.
// Return one randomly if there is more than one point.
//
// Check steps 3-4 in Algorithm 3 in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
size_t NSGAIII<TypeIndiv,TypePop>::FindNicheReferencePoint(const vector<ReferencePoint> &rps)
{
	// find the minimal cluster size
	size_t min_size = numeric_limits<size_t>::max();
	for (size_t r=0; r<rps.size(); r+=1)
	{
		min_size = std::min(min_size, rps[r].MemberSize());
	}

	// find the reference points with the minimal cluster size Jmin
	vector<size_t> min_rps;
	for (size_t r=0; r<rps.size(); r+=1)
	{
		if (rps[r].MemberSize() == min_size)
		{
			min_rps.push_back(r);
		}
	}
	
	// return a random reference point (j-bar)
	int rnd=Global::msp_global->getRandInt(0,min_rps.size());
	return min_rps[rnd];
}

// ----------------------------------------------------------------------
// SelectClusterMember():
//
// Select a potential member (an individual in the front Fl) and associate
// it with the reference point.
//
// Check the last two paragraphs in Section IV-E in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
int NSGAIII<TypeIndiv,TypePop>::SelectClusterMember(const ReferencePoint &rp)
{
	int chosen = -1;
	if (rp.HasPotentialMember())
	{
		if (rp.MemberSize() == 0) // currently has no member
		{
			chosen =  rp.FindClosestMember();
		}
		else
		{
			chosen =  rp.RandomMember();
		}
	}
	return chosen;
} 

// ----------------------------------------------------------------------
// EnvironmentalSelection():
//
// Check Algorithms 1-4 in the original paper.
// ----------------------------------------------------------------------
template<typename TypeIndiv,typename TypePop>
void NSGAIII<TypeIndiv,TypePop>::EnvironmentalSelection()
{
	vector<ReferencePoint> rps=mv_rps;
	vector<vector<int> > fronts;
	int rank=0;
	int count=0;
	int size=m_offspring.getPopSize();
	// ---------- Step 4 in Algorithm 1: non-dominated sorting ----------
	m_offspring.rank();
	while(1)
	{
		vector<int> temp;
		for(int i=0;i<size;i++)
		{
			if(m_offspring[i]->rank()==rank)
			{
				temp.push_back(i);
				++count;
			}
		}
		fronts.push_back(temp);
		if(count==size) break;
		++rank;
	}
	
	// ---------- Steps 5-7 in Algorithm 1 ----------
	vector<size_t> considered; // St
	int last = 0, next_size = 0;
	while (next_size < m_parent.getPopSize())
	{
		next_size += fronts[last].size();
		last += 1; 
	}
	fronts.erase(fronts.begin()+last, fronts.end()); // remove useless individuals

	count=0;
	for (size_t t=0; t<fronts.size()-1; t+=1)
		for (size_t i=0; i<fronts[t].size(); i+=1)
			*m_parent[count++]=*m_offspring[fronts[t][i]];

	// ---------- Steps 9-10 in Algorithm 1 ----------
	if (count == m_parent.getPopSize()) return;


	// ---------- Step 14 / Algorithm 2 ----------
	vector<double> ideal_point = TranslateObjectives(fronts);

	vector<size_t> extreme_points;
	FindExtremePoints(&extreme_points, fronts);

	vector<double> intercepts;
	ConstructHyperplane(&intercepts, extreme_points);

	NormalizeObjectives(fronts, intercepts, ideal_point);

	// ---------- Step 15 / Algorithm 3, Step 16 ----------
	Associate(&rps, mvv_offConvObj , fronts);

	// ---------- Step 17 / Algorithm 4 ----------
	while (count < m_parent.getPopSize())
	{
		size_t min_rp = FindNicheReferencePoint(rps);

		int chosen = SelectClusterMember(rps[min_rp]);
		if (chosen < 0) // no potential member in Fl, disregard this reference point
		{
			rps.erase(rps.begin()+min_rp); 
		}
		else
		{
			rps[min_rp].AddMember();
			rps[min_rp].RemovePotentialMember(chosen);
			*m_parent[count++]=*m_offspring[chosen];
		}
	}

}
// ----------------------------------------------------------------------

#endif