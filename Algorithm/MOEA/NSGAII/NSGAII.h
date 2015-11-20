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

*  See the details of NSGA2 in the following paper
*  A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II 
*  Kalyanmoy Deb, Associate Member, IEEE, Amrit Pratap, Sameer Agarwal, and T. Meyarivan
*  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 6, NO. 2, APRIL 2002
*************************************************************************/
// Created: 7 Jan 2015
// Last modified:

#ifndef NSGAII_H
#define NSGAII_H

#include "../../DE/DEPopulation.h"
#include "../../../Measure/mMultiObj.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

template<typename TypeIndiv,typename TypePop>   
class NSGAII : public Algorithm
{
public:
	NSGAII(ParamMap &v);
	~NSGAII();
	ReturnFlag run_();

protected:
	void eval_dens();
	virtual void evolve_mo()=0;
	int tour_selection();
	TypePop m_parent, m_offspring;
};

template<typename TypeIndiv,typename TypePop>
NSGAII<TypeIndiv,TypePop>::NSGAII(ParamMap &v): m_parent(v[param_popSize]),m_offspring(2*(v[param_popSize]),false),Algorithm(-1,string())
{
	for(int i=0;i<m_parent.getPopSize();i++)
		m_parent[i]->rank()=0;
	m_parent.clearBestArchive();
	m_offspring.clearBestArchive();
}

template<typename TypeIndiv,typename TypePop>
NSGAII<TypeIndiv,TypePop>::~NSGAII()
{

}

template<typename TypeIndiv,typename TypePop>
int NSGAII<TypeIndiv,TypePop>::tour_selection()
{
	int p1 = Global::msp_global->getRandInt(0,m_parent.getPopSize());
	int p2 = Global::msp_global->getRandInt(0,m_parent.getPopSize());

	if(m_parent[p1]->rank()<m_parent[p2]->rank())
		return p1;
	else
		return p2;
}

template<typename TypeIndiv,typename TypePop>
void NSGAII<TypeIndiv,TypePop>::eval_dens()
{
	int numobj=Global::msp_global->mp_problem->getNumObj();
	int pops=0;  //indicate parent population size be 0
	int size = m_offspring.getPopSize();
	int rank = 0;
	while(1){
		int count = 0;
        for(int i=0; i<size; i++)
			if(m_offspring[i]->rank()==rank)
				count++;

		int size2 = pops + count;
		if(size2>m_parent.getPopSize()) {
			break;
		}

        for(int i=0; i<size; i++)
  	        if(m_offspring[i]->rank()==rank)
			{
				*m_parent[pops]=*m_offspring[i];
				++pops;
			}

		rank++;
		if(pops>=m_parent.getPopSize()) break;
	}

	if(pops<m_parent.getPopSize()){
		vector<int> list;
		// save the individuals in the overflowed front
        for(int i=0; i<size; i++)
  	        if(m_offspring[i]->rank()==rank)
				list.push_back(i);
		int s2 = list.size();
		vector<double> density(s2);
		vector<double> obj(s2);
		vector<int> idx(s2);
		vector<int> idd(s2);
	
		for(int i=0; i<s2; i++){
			idx[i]     = i;
			density[i] = 0;
		}

		for(int j=0; j<numobj; j++){		    			
			for(int i=0; i<s2; i++){
			    idd[i] = i;
				obj[i] = m_offspring[list[i]]->data().m_obj[j];
			}
			//gMinfastsort(obj,idd,s2,s2);
			gQuickSort(obj,s2,idd,true,0,s2-1,s2);
            density[idd[0]]    += -1.0e+30;
            density[idd[s2-1]] += -1.0e+30;
			for(int k=1; k<s2-1; k++)
				density[idd[k]]+= -(obj[idd[k]] - obj[idd[k-1]] + obj[idd[k+1]] - obj[idd[k]]);
		}
		idd.clear();
		obj.clear();

		int s3 = m_parent.getPopSize() - pops;

		//gMinfastsort(density,idx,s2,s3);
		gQuickSort(density,s2,idx,true,0,s2-1,s3);
		for(int i=0; i<s3; i++)
		{
			*m_parent[pops]=*m_offspring[list[idx[i]]];
			++pops;
		}

		density.clear();
		idx.clear();
		list.clear();
	}
}

template<typename TypeIndiv,typename TypePop>
ReturnFlag NSGAII<TypeIndiv,TypePop>::run_()
{
 	
	#ifdef OFEC_CONSOLE
	if(mMultiObj::getMultiObj()&&Global::msp_global->mp_problem->isGlobalOptKnown())
		mMultiObj::getMultiObj()->recordDistance<TypeIndiv>(Global::msp_global.get(),Global::msp_global->m_runId,m_parent.getPop());
#endif

	// evolution
	while(!this->ifTerminating())
	{
	//	cout<<"Run "<<Global::msp_global->m_runId<<" is running "<<Global::msp_global->mp_problem->getEvaluations()<<endl;
	
		evolve_mo();
		m_offspring.rank();
		eval_dens();

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
		mMultiObj::getMultiObj()->record(Global::msp_global->m_runId,i,m_parent[i]->data().m_obj,m_parent[i]->data().m_x);
	}
#endif

//	cout<<"Run "<<Global::msp_global->m_runId<<" is terminated"<<endl;
	return Return_Normal;
		
}


#endif //NSGAII_H