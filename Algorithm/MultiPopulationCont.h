/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
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
// Created: 21 September 2011
// Last modified: 12 Dec. 2014
#ifndef MULTI_POPULATION_CONT_H
#define MULTI_POPULATION_CONT_H

#include "PopulationCont.h"
#include "MultiPopulation.h"

template<typename TypePop>
class MultiPopulationCont:public MultiPopulation<TypePop>{
protected:
	double m_overlapDegree;
public:
	virtual ~MultiPopulationCont(){}
	MultiPopulationCont():MultiPopulation<TypePop>(){}
	MultiPopulationCont(int size,double degree=0):MultiPopulation<TypePop>(size),m_overlapDegree(degree){}
	MultiPopulationCont(const int num,const int subSize,bool mode, const int runId):MultiPopulation<TypePop>(num,subSize,mode,runId){}
	void printPops(ofstream &out,Problem *pro=0);
	double getAvgInitialRadius();
	double getAvgCurRadius();
	double getAvgRadiusQaulity(Problem *pro);
	bool isAvgRadiusBelow(double percent, double threshold);
	double getAvgDistanceBetwPop();
	virtual int removeOverlapping();
	void setOverlapDgre(double degree);
	void handleReturnFlagAll(ReturnFlag f);
	void measureMultiPop();
};

template<typename TypePop>
void MultiPopulationCont<TypePop>::printPops(ofstream &out,Problem *pro){	
	//print peaks
	dynamic_cast<DynamicContinuous*> (pro)->printPeaks(out);
	//out<<endl;

	/// print centers of each pop in gnuplot format for 2-D problems
	for(unsigned int i=0;i<this->m_subPop.size();i++){
		Solution<CodeVReal> &x=this->m_subPop[i]->m_center;
		out<<x.data()[0]<<" "<<x.data()[1]<<" "<<this->m_subPop[i]->m_curRadius<<endl;
		//out<<"plot [0:2*pi] "<<this->m_subPop[i]->m_curRadius<<"*sin(t)+"<<x.getGene<double>(0)<<","<<this->m_subPop[i]->m_curRadius<<"*cos(t)+"<<x.getGene<double>(1)<<" with lines lt 1"<<endl;
		//out<<"plot [0:2*pi] "<<this->m_subPop[i]->m_initialRadius<<"*sin(t)+"<<x.getGene<double>(0)<<","<<this->m_subPop[i]->m_initialRadius<<"*cos(t)+"<<x.getGene<double>(1)<<" with lines lt 2"<<endl;
	}
	out<<"#########line"<<dynamic_cast<DynamicContinuous*> (pro)->getNumberofPeak()+this->m_subPop.size()<<"################################"<<endl;
	MultiPopulation<TypePop>::printPops(out,pro);	
}


template<typename TypePop>
double MultiPopulationCont<TypePop>::getAvgInitialRadius(){
	if(this->m_subPop.size()==0) return 0;
	double r=0;
	for(unsigned int j=0;j<this->m_subPop.size();j++){
		r+=this->m_subPop[j]->m_initialRadius;
	}
	r/=this->m_subPop.size();
	return r;
}

template<typename TypePop>
double MultiPopulationCont<TypePop>::getAvgCurRadius(){
	if(this->m_subPop.size()==0) return 0;
	double r=0;
	for(unsigned int j=0;j<this->m_subPop.size();j++){
		r+=this->m_subPop[j]->m_curRadius;
	}
	r/=this->m_subPop.size();
	return r;
}

template<typename TypePop>
double MultiPopulationCont<TypePop>::getAvgRadiusQaulity(Problem *pro){
	if(this->m_subPop.size()==0) return 0;
	for(unsigned j=0;j<this->m_subPop.size();j++){
		if(this->m_subPop[j]->m_popsize==0) continue;
		double minDis;int nearest=-1,numPeaksIn=0;
		// get the neaest peak which allocates within its radius
		for(int k=0;k<dynamic_cast<DynamicProblem*>(pro)->getNumberofPeak();k++ ){
			if(!dynamic_cast<DynamicContinuous*>(pro)->isVisable(k)) continue;
			const double *p=dynamic_cast<DynamicContinuous*>(pro)->getPeak(k);
			CodeVReal peak(p,p+GET_NUM_DIM);
			double dis=this->m_subPop[j]->getNearestBest2Peak(peak)->getDistance(peak);
			if(dis<this->m_subPop[j]->m_initialRadius){
				numPeaksIn++;
				if(nearest==-1){
					minDis=dis;
					nearest=k;
				}else{
					if(minDis>dis){
						minDis=dis;
						nearest=k;
					}
				}
			}
		}
		if(nearest!=-1) this->m_subPop[j]->updateRadiusQaulity(nearest,numPeaksIn,this->m_subPop[j]->m_initialRadius);
	}
	double rq=0.; int count=0;
	for(unsigned int j=0;j<this->m_subPop.size();j++){

		if(this->m_subPop[j]->m_radiusQaulity!=0){
			rq+=this->m_subPop[j]->m_radiusQaulity;
			count++;
		}
	}
	if(count) rq/=count;
	else rq=0;
	return rq;
}

template<typename TypePop>
bool MultiPopulationCont<TypePop>::isAvgRadiusBelow(double percent, double threshold){
	int count=0;
	for(int j=0;j<this->m_subPop.size();j++){
		if(this->m_subPop[j]->m_curRadius<threshold) count++;
	}
	if(count*1./this->m_subPop.size()>=percent) return true;
	else return false;
}

template<typename TypePop>
double MultiPopulationCont<TypePop>::getAvgDistanceBetwPop(){
	if(this->m_subPop.size()==0) return 0;
	double r=0;
	for(unsigned int j=0;j<this->m_subPop.size();j++){
		for(unsigned int i=j+1;i<this->m_subPop.size();i++)
		r+=this->m_subPop[j]->m_center.getDistance(this->m_subPop[i]->m_center);
	}
	r=r/((this->m_subPop.size()+1)*this->m_subPop.size()/2.);
	return r;
}

template<typename TypePop>
int MultiPopulationCont<TypePop>::removeOverlapping(){
	for(unsigned i=0;i<this->m_subPop.size();i++){
		if(this->m_subPop[i]->m_popsize==0) continue;
		for(unsigned j=i+1;j<this->m_subPop.size();j++){	
			if(this->m_subPop[j]->m_popsize==0) continue;
			double dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_center);
			if(dist<this->m_subPop[i]->m_initialRadius||dist<this->m_subPop[j]->m_initialRadius){
				int c1=0,c2=0;
				for(int k=0;k<this->m_subPop[j]->m_popsize;k++){
					dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c1++;
				}
				for(int k=0;k<this->m_subPop[i]->m_popsize;k++){
					dist=this->m_subPop[j]->m_center.getDistance(this->m_subPop[i]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c2++;
				}
				if(c1>this->m_subPop[j]->m_popsize*m_overlapDegree&&c2>this->m_subPop[i]->m_popsize*m_overlapDegree){
					int idx=-1;
					if(*this->m_subPop[i]>(*this->m_subPop[j])){	
						this->m_subPop[i]->add(*this->m_subPop[j]);
						this->deletePopulation(j);
						idx=j;
					}else{
						this->m_subPop[j]->add(*this->m_subPop[i]);
						this->deletePopulation(i);
						idx=i;
					}
					return idx;
				}
			}
		}
	}
	return -1;
}

template<typename TypePop>
void MultiPopulationCont<TypePop>::setOverlapDgre(double degree){
	m_overlapDegree=degree;
}

template<typename TypePop>
void MultiPopulationCont<TypePop>::handleReturnFlagAll(ReturnFlag f){
	MultiPopulation<TypePop>::handleReturnFlagAll(f);
	if(mMultiModal::getPopInfor()){
	if(f==Return_ChangeNextEval){
		measureMultiPop();	
	}
	}
}

template<typename TypePop>
void MultiPopulationCont<TypePop>::measureMultiPop(){
	if(Global::msp_global->mp_algorithm->ifTerminated()) return;
	if(mMultiModal::getPopInfor()){
	
	if(IS_PROBLEM_NAME(Global::ms_curProId,"DYN_CONT_MovingPeak")){
		int peaksf=CAST_PROBLEM_DYN_CONT->getPeaksFound();
		mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
			Global::msp_global->m_totalNumIndis,this->getNumPops()+this->m_convergedPops,peaksf,\
			CAST_PROBLEM_DYN_CONT->getNumofVisablePeaks(),getAvgInitialRadius(),getAvgCurRadius(),0,0,\
			CAST_PROBLEM_DYN_CONT->getPeaksTracedQaulity(),getAvgRadiusQaulity(Global::msp_global->mp_problem.get()),\
																	CAST_PROBLEM_DYN_CONT->isGOptTracked());
	}else if(IS_PROBLEM_NAME(Global::ms_curProId,"DYN_CONT_RotationDBG")){
		int peaksf=CAST_PROBLEM_DYN_CONT->getPeaksFound();
		mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
			Global::msp_global->m_totalNumIndis,this->getNumPops()+this->m_convergedPops,peaksf,\
			CAST_PROBLEM_DYN_CONT->getNumofVisablePeaks(),getAvgInitialRadius(),getAvgCurRadius(),0,0,\
			0,0,CAST_PROBLEM_DYN_CONT->isGOptTracked());
	}/*else if(IS_PROBLEM_NAME(Global::ms_curProId, "FUN_FreePeak_D_OnePeak")){
		FFreePeak_D_OnePeak*pro = dynamic_cast<FFreePeak_D_OnePeak*>(Global::msp_global->mp_problem.get());
		int peaksf = pro->getPT();
		mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(), \
			Global::msp_global->m_totalNumIndis, this->getNumPops() + this->m_convergedPops, peaksf, \
			pro->getNumPeak(), getAvgInitialRadius(), getAvgCurRadius(), 0, 0, \
			0, 0, pro->isAllGOptTraced());
	}*/
	else if (Global::msp_global->mp_problem->m_name.find("FUN_") != string::npos){
			int peaksf=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
		    mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
			Global::msp_global->m_totalNumIndis,this->getNumPops()+this->m_convergedPops,peaksf,\
			CAST_PROBLEM_CONT->getGOpt().getNumOpt(),getAvgInitialRadius(),getAvgCurRadius(),0,0,\
			0,0,CAST_PROBLEM_CONT->getGOpt().isAllFound());
	}

	}
}
#endif