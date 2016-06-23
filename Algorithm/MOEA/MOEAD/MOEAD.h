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
// Created: 30 December 2014
// Last modified:

#ifndef MOEAD_H
#define MOEAD_H

#include "../../DE/DEPopulation.h"
#include "../../../Measure/mMultiObj.h"
#include "../../../Utility/myVector.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

template<typename TypeIndiv,typename TypePop>
class MOEAD: public TypePop
{
	enum DecomFun{_TCHE1,_TCHE2,_NBI1,_NBI2,_NBI3};
public:
	MOEAD(ParamMap &v);
	~MOEAD();
	virtual void evolve_mo()=0;
	ReturnFlag run_();
	void setDecomFun(DecomFun f);
protected:
	void init_uniformweight();
	void init_neighbourhood();
	void update_reference(TypeIndiv &sol);
	void update_problem(TypeIndiv &sol, int id, int type);
	double fitnessfunction(vector<double> &obj,int k);
	void matingselection(vector<int> &list, int cid, int size, int type);
	vector<double> mv_idealPoint; //the best value in every dimension
	vector<TypeIndiv> mv_solArr; //corresponding TypeIndiv of the best value in every dimension
	vector<MyVector > mv_namda;
	vector<vector<int> > mvv_neigh;
	int m_unit;
	int m_limit;
	int m_niche; //number of neighbours
	double m_realb;     // probability of selecting mating parents from neighborhood
	DecomFun m_decomFunction;
};


template<typename TypeIndiv,typename TypePop>
MOEAD<TypeIndiv,TypePop>::MOEAD(ParamMap &v):TypePop(v[param_popSize]),m_unit(33),m_niche(20),m_realb(0.9),m_limit(2),m_decomFunction(_TCHE1)
{
	int numObj=Global::msp_global->mp_problem->getNumObj();
	mv_solArr.resize(numObj);
	if(Global::msp_global->mp_problem->getOptType()==MIN_OPT)
		mv_idealPoint.resize(numObj,1.0e+30);
	else
		mv_idealPoint.resize(numObj,-1.0e+30);
	init_uniformweight();
	init_neighbourhood();
	if(this->m_popsize!=this->m_pop.size())
	{
		this->m_popsize=mv_namda.size();
		this->m_pop.clear();
		this->m_pop.resize(this->m_popsize);
		for(auto &i:this->m_pop) i=move(unique_ptr<TypeIndiv>(new TypeIndiv()));
		Global::msp_global->mp_problem->resetEvaluations();
		this->initialize(false,true,true);
	}
	for(int i=0;i<this->m_popsize;i++)
		update_reference(*(this->m_pop[i]));
}

template<typename TypeIndiv,typename TypePop>
MOEAD<TypeIndiv,TypePop>::~MOEAD()
{
	mv_idealPoint.clear();
	mv_solArr.clear();
	mv_namda.clear();
	mvv_neigh.clear();
}


template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::init_uniformweight()
{
	
	if(Global::msp_global->mp_problem->getNumObj()==2)
	{
		mv_namda.resize(this->m_popsize);
		for(int n=0; n<this->m_popsize; n++)
		{
			double a = 1.0*n/(this->m_popsize - 1);
			mv_namda[n].push_back(a);
			mv_namda[n].push_back(1-a);
		}
	}
	else
	{
		int n=0;
		for(int i=0; i<=m_unit; i++)
		{
			for(int j=0; j<=m_unit; j++)
			{
				if(i+j<=m_unit)
				{
					vector<int> arr;
					arr.push_back(i);
					arr.push_back(j);
					arr.push_back(m_unit-i-j);
					mv_namda.push_back(vector<double>(0));
					for(int k=0; k<arr.size(); k++)
						mv_namda[n].push_back(1.0*arr[k]/m_unit);
					n++;
					arr.clear();
				}
			}
		}
		this->m_popsize=n;
	}
}

template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::init_neighbourhood()
{
	int pops=mv_namda.size();
	mvv_neigh.resize(pops);
	vector<double> dis(pops);
	vector<int> index(pops);
	for(int i=0; i<pops; i++)
	{
		// calculate the distances based on weight vectors
		for(int j=0; j<pops; j++)
		{
			dis[j] = mv_namda[i].getDis(mv_namda[j]);
			index[j] = j;
		}
	
		// find 'niche' nearest neighboring subproblems
		//gMinfastsort(dis,index,pops,m_niche);
		gQuickSort(dis,pops,index,true,0,pops-1,m_niche);
		for(int k=0; k<m_niche; k++)
			mvv_neigh[i].push_back(index[k]);
	}
	dis.clear();
	index.clear();
}

template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::update_reference(TypeIndiv &sol)
{
	//sol: child TypeIndiv
	int numObj=Global::msp_global->mp_problem->getNumObj();
	if(Global::msp_global->mp_problem->getOptType()==MIN_OPT)
	{
		for(int n=0; n<numObj; n++)
		{
			if(sol.data().m_obj[n]<mv_idealPoint[n])
			{
				mv_idealPoint[n] = sol.data().m_obj[n];
				mv_solArr[n]    = sol;
			}
		}
	}
	else
	{
		for(int n=0; n<numObj; n++)
		{
			if(sol.data().m_obj[n]>mv_idealPoint[n])
			{
				mv_idealPoint[n] = sol.data().m_obj[n];
				mv_solArr[n]    = sol;
			}
		}
	}
}


template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::update_problem(TypeIndiv &sol, int id, int type)
{
	// sol: child TypeIndiv
	// id:   the id of current subproblem
	// type: update TypeIndivs in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if(type==1)	size = m_niche;
	else        size = this->m_popsize;
	vector<int> perm(size);
	Global::msp_global->initializeRandomArray<vector<int> >(perm,size);
    for(int i=0; i<size; i++)
	{
		int k;
		if(type==1) k = mvv_neigh[id][perm[i]];
		else        k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		f1 = fitnessfunction(this->m_pop[k]->data().m_obj, k);
		f2 = fitnessfunction(sol.data().m_obj, k);
		if(f2<f1)
		{
			*this->m_pop[k] = sol;
			time++;
		}
		// the maximal number of TypeIndivs updated is not allowed to exceed 'limit'
		if(time>=m_limit)
			return;
	}
	perm.clear();
}

template<typename TypeIndiv,typename TypePop>
double MOEAD<TypeIndiv,TypePop>::fitnessfunction(vector<double> &obj,int k)
{
	// Chebycheff Scalarizing Function
	double fitness = 0;
	int numObj=obj.size();
	if(m_decomFunction==_TCHE1)
	{
		double max_fun = -1.0e+30;
		for(int n=0; n<numObj; n++)
		{
			//double diff = fabs(y_obj[n] - idealpoint[n] + scale[n]);
			//double diff = fabs(y_obj[n] - idealpoint[n] + 0.05);
		    double diff = fabs(obj[n] - mv_idealPoint[n]);
			//double diff = fabs(y_obj[n] - 0);
			double feval;
			if(mv_namda[k][n]==0) 
				feval = 0.0001*diff;
			else
			    feval = diff*mv_namda[k][n];
			if(feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}

	if(m_decomFunction==_TCHE2)
	{
		// reference point in the CHIM
		vector<int> scale(numObj);
		throw myException("Please initialize the scale @MOEAD<TypeIndiv,TypePop>::fitnessfuction");
		double max_fun = -1.0e+30;
		for(int n=0; n<numObj; n++)
		{
			double diff = (obj[n] - mv_idealPoint[n])/scale[n];  //note: the scale is not initialized, there has no knowledge
			double feval;
			if(mv_namda[k][n]==0) 
				feval = 0.0001*diff;
			else
			    feval = diff*mv_namda[k][n];
			if(feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}
	


	// CHIM + Tchebycheff
	// CHIM is not available in 3 objectives
	if(m_decomFunction==_NBI1){

		// quasi normal direction
		MyVector norm;
		for(int i=0; i<numObj; i++)
		{
		   norm.push_back(0.0);
		   for(int j=0; j<numObj; j++){
			   norm[i]+= -mv_solArr[j].data().m_obj[i];
		   }
		}

		// normalization
		norm.normalize();

		// reference point in the CHIM
		vector <double> base;
		for(int i=0; i<numObj; i++)
		{
			double tp2 = 0;
			for(int j=0; j<numObj; j++)
				tp2+= mv_solArr[j].data().m_obj[i]*mv_namda[k][j];
			base.push_back(tp2);
		}

		// Tchebycheff function
		double max_fun = -1.0e+30;
		for(int n=0; n<numObj; n++)
		{
			double diff  = obj[n] - base[n];
			double feval = -diff*norm[n];
			if(feval>max_fun) max_fun = feval;

		}	
		fitness = max_fun;
	}

	//* Boundary intersection approach
	//* reference point is chosen as the ideal point
	//* the direction is independent of CHIM
	if(m_decomFunction==_NBI2)
	{

		mv_namda[k].normalize();

	    // penalty method 
	    // temporary vectors NBI method
		MyVector realA(numObj);
		MyVector realB(numObj);

		// difference beween current point and reference point
		for(int n=0; n<numObj; n++)
			realA[n] = (obj[n] - mv_idealPoint[n]);

		// distance along the search direction norm
		double d1 = fabs(realA*mv_namda[k]);

		// distance to the search direction norm
		for(int n=0; n<numObj; n++)
			realB[n] = (obj[n] - (mv_idealPoint[n] + d1*mv_namda[k][n]));
		double d2 = realB.length();

		fitness =  (d1 + 5*d2);

		//t2 = clock();
	    //total_sec+=(t2 - t1);
	}

	// NBI method
	if(m_decomFunction==_NBI3){

		// quasi normal direction
		MyVector norm;
		for(int i=0; i<numObj; i++)
		{
		   norm.push_back(0.0);
		   for(int j=0; j<numObj; j++){
			   norm[i]+= -mv_solArr[j].data().m_obj[i];
		   }
		}

		// normalization
		norm.normalize();


		// reference point in the CHIM
		vector <double> base;
		for(int i=0; i<numObj; i++)
		{
			double tp2 = 0;
			for(int j=0; j<numObj; j++)
				tp2+= mv_solArr[j].data().m_obj[i]*mv_namda[k][j];
			base.push_back(tp2);
		}

	    // penalty method 
	    // temporary vectors NBI method
		MyVector realA;
		MyVector realB;

		// difference beween current point and reference point
		for(int n=0; n<numObj; n++)
			realA.push_back(obj[n] - base[n]);

		// distance along the search direction norm
		double d1 = realA*norm;

		// distance to the search direction norm
		for(int n=0; n<numObj; n++)
			realB.push_back(obj[n] - (base[n] + d1*norm[n]));
		double d2 = realB.length();

		fitness =  -d1 + 2*d2;
	}
	return fitness;
}

template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::matingselection(vector<int> &list, int cid, int size, int type)
{
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int r, p;
    while(list.size()<size)
	{
		if(type==1){
			r = Global::msp_global->getRandInt(0,m_niche);
			p = mvv_neigh[cid][r];
		}
		else
			p = Global::msp_global->getRandInt(0,this->m_popsize);

		bool flag = true;
		for(int i=0; i<list.size(); i++)
		{
			if(list[i]==p) // p is in the list
			{ 
				flag = false;
				break;
			}
		}

		if(flag) list.push_back(p);
	}
}


template<typename TypeIndiv,typename TypePop>
ReturnFlag MOEAD<TypeIndiv,TypePop>::run_()
{
	#ifdef OFEC_CONSOLE
	if(mMultiObj::getMultiObj()&&Global::msp_global->mp_problem->isGlobalOptKnown())
		mMultiObj::getMultiObj()->recordDistance<TypeIndiv>(Global::msp_global.get(),Global::msp_global->m_runId,this->m_pop);
	#endif

	// evolution
	while(!this->ifTerminating())
	{
		//cout << "Run " << Global::msp_global->m_runId << "  " << Global::msp_global->mp_problem->getEvaluations() << " " << mMultiObj::getMultiObj()->getCurDis2PF(Global::msp_global->m_runId) << endl;


	
		evolve_mo();
#ifdef OFEC_DEMON
		vector<Algorithm*> vp;
		this->rank();
		vp.push_back(this);
		msp_buffer->updateBuffer_(&vp);
#endif
		#ifdef OFEC_CONSOLE
		if(mMultiObj::getMultiObj()&&Global::msp_global->mp_problem->isGlobalOptKnown())
			mMultiObj::getMultiObj()->recordDistance<TypeIndiv>(Global::msp_global.get(),Global::msp_global->m_runId,this->m_pop);
#endif

	}
	#ifdef OFEC_CONSOLE
	if(mMultiObj::getMultiObj()){
		mMultiObj::getMultiObj()->reInitialize(Global::msp_global.get(),this->m_popsize);
		for(int i=0;i<this->m_popsize;i++)
			mMultiObj::getMultiObj()->record(Global::msp_global.get(),i,this->m_pop[i]->data().m_obj,this->m_pop[i]->data().m_x);
	}
#endif
//	cout<<"Run "<<Global::msp_global->m_runId<<" is terminated"<<endl;
	return Return_Normal;
}

template<typename TypeIndiv,typename TypePop>
void MOEAD<TypeIndiv,TypePop>::setDecomFun(DecomFun f){
	m_decomFunction=f;
}
#endif //MOEAD_H
