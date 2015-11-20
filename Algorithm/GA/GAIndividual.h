/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
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
// Created: 12 Nov 2011
// Last modified:

#ifndef GAINDIVIDUAL_H
#define GAINDIVIDUAL_H

#include "../Individual.h"
#include "../../Global/boundary.h"
#include "../../Problem/ContinuousProblem.h"

template<typename, typename> class GAPopulation;
class SOSSubPop;

template<typename ED>
class GAIndividual : public Individual<ED>{
	template<typename, typename> friend class GAPopulation;
	friend class SOSSubPop;
 public:
        GAIndividual();
		virtual void xover(GAIndividual<ED> & mate);

        virtual ~GAIndividual();
        GAIndividual(const GAIndividual<ED> & );
		GAIndividual<ED>& operator=(const GAIndividual<ED> & rhs);
		virtual void mutate();
		
        void setDistribution(double dis);
        void setMutationProbability(double p);
        void setXoverProbability(double p);
        void increaseDimension();
        void decreaseDimension();
		void setMutationStrategy(GAMutationStrategy s);
		void setMutationStep(double step);
		void setXoverStrategy(GAXoverStrategy s);
		void mutateReal();
protected:
        Solution<ED> m_backup;            // backup individual for p_self
        double m_distribution;          // distribution for mutation for real-coded encoding
        double m_mutationProb;          // mutation probability
        double m_xoverProb;             //crossover probability
		GAMutationStrategy m_mutationStrategy;
		GAXoverStrategy m_xoverStrategy;
		double m_mutationStep;			// mutation step in normal mutation strategy     
    protected:
		void polynomialMutate(const DomainReal &domain,int site);
		void normalMutate(int);      
        void realXover(GAIndividual<ED> & mate);
        double getDelta(double u, double delta_l, double delta_u);
		void singlePointXover(GAIndividual<ED> & mate);
		template<typename DT>
		void mutateOneGene(const DT &, int );
};

template<typename ED>
GAIndividual<ED>::GAIndividual():Individual<ED>(),m_distribution(100),m_mutationProb(0.01),m_xoverProb(0.8),m_mutationStrategy(MUTAT_NORMAL),m_mutationStep(0.1),\
	m_xoverStrategy(XOVER_ARITHMETIC){
    //ctor

}

template<typename ED> 
GAIndividual<ED>::~GAIndividual()
{
    //dtor

}

template<typename ED> 
GAIndividual<ED>::GAIndividual(const GAIndividual<ED> & rhs):Individual<ED>(rhs){

     m_distribution=rhs.m_distribution;
     m_mutationProb=rhs.m_mutationProb;
     m_xoverProb=rhs.m_xoverProb;

}

template<typename ED>
 GAIndividual<ED>&  GAIndividual<ED>::operator=(const GAIndividual<ED> & rhs){
     if(this==&rhs) return *this;
     Individual<ED>::operator=(rhs);

     m_distribution=rhs.m_distribution;
     m_mutationProb=rhs.m_mutationProb;
     m_xoverProb=rhs.m_xoverProb;
    
     return *this;

 }

 template<typename ED>
 void  GAIndividual<ED>::setMutationStrategy(GAMutationStrategy s){
	 m_mutationStrategy=s;

 }

 template<typename ED>
 void  GAIndividual<ED>::setXoverStrategy(GAXoverStrategy s){
	 m_xoverStrategy=s;
 }

 template<typename ED>
 void  GAIndividual<ED>::setMutationStep(double step){
	 m_mutationStep=step;

 }
 
template<typename ED> template<typename DT>
void GAIndividual<ED>::mutateOneGene(const DT  &domain, int site){
	 
	 if(typeid(DT)==typeid(DomainReal)){
			if(m_mutationStrategy==MUTAT_POLYNOMIAL){
				polynomialMutate(domain,site);
			}else if(m_mutationStrategy==MUTAT_NORMAL){
				normalMutate(site);
			}
	 }else{
		if(m_mutationStrategy==MUTAT_COMB){	
			//note: the value may not be mutated
			this->data()[site]=Global::msp_global->getRandInt(domain.m_lower,domain.m_upper);
		}
	 }
	 
 }
 
template<typename ED>
void  GAIndividual<ED>::mutate(){
    //operation performed must be after crossover
	mutateReal();
}

template<typename ED>
void GAIndividual<ED>::mutateReal(){
	bool flag=false;
	int numDim=GET_NUM_DIM;
	for(int i=0;i<numDim;i++){
		DomainReal &bound=CAST_PROBLEM_CONT->getSearchRange().getDomain(i);		
		if(Global::msp_global->mp_uniformAlg->Next()<=m_mutationProb){
			flag=true;
			mutateOneGene(bound,i);
		}	
	}
	if(flag) this->self().validate();
}

template<typename ED>
void GAIndividual<ED>::normalMutate(int site){
	this->data()[site]=this->data()[site]+m_mutationStep*Global::msp_global->mp_normalAlg->Next();
}

template<typename ED>
void  GAIndividual<ED>::polynomialMutate(const DomainReal &domain,int site){
//===================================================================
//Mutation Using polynomial probability distribution. Picks up a random
//site and generates a random number u between -1 to 1, ( or between
//minu to maxu in case of rigid boudaries) and calls the routine
//get_delta() to calculate the actual shift of the value.
//
//This is from http://www.iitk.ac.in/kangal/codes.shtml
//:Single-objective GA code in C (for Windows and Linux).
//    ====================================================================
	double distance1,x,delta_l,delta_u,delta,u,upper,lower;

	if(CAST_PROBLEM_CONT->getBoundaryFlag(site)){
			
		lower=domain.m_lower; upper=domain.m_upper;
		x = this->data()[site];
			
		distance1 = lower - x;
		delta_l = distance1/(upper - lower);
		if (delta_l < -1.0)  delta_l = -1.0;

		distance1 = upper - x;
		delta_u = distance1/(upper - lower);
		if (delta_u > 1.0)   delta_u = 1.0;

		if (-1.0*delta_l < delta_u) delta_u = -1.0 * delta_l;
		else delta_l = -1.0 * delta_u;
	}else{
		delta_l = -1.0;
		delta_u =  1.0;
	}

	u = Global::msp_global->mp_uniformAlg->Next();
	// calculation of actual delta value 
	delta = getDelta(u, delta_l, delta_u) * (upper - lower);

	this->data()[site] =this->data()[site]+ delta;
			
}

template<typename ED>
void GAIndividual<ED>::xover(GAIndividual<ED> & mate){
//====================================================================
//CROSS - OVER  USING strategy of uniform 50% variables
//  For one variable problem (each gene contains one variable), it is crossed over as usual.
//  For multivariables, each variable is crossed over with a probability
//  of 50 % , each time generating a new random beta.
//====================================================================
	switch (m_xoverStrategy)
	{
		case XOVER_ARITHMETIC:
			 realXover(mate);
			break;
		case XOVER_SINGLEPOINT:
			singlePointXover(mate);
			break;
		default:
			throw myException("crossover not defined!@template<typename ED> GAIndividual<ED>::xover");
	}
  
}

template<typename ED>
double GAIndividual<ED>::getDelta(double u, double delta_l, double delta_u)
//==================================================================
//For given u value such that   -1 <= u <= 1, this routine returns a
//value of delta from -1 to 1. Exact value of delta depends on specified
//n_distribution.
//
//This is from http://www.iitk.ac.in/kangal/codes.shtml
//:Single-objective GA code in C (for Windows and Linux).
//====================================================================
{
  double delta, aa;

  if (u >= 1.0-1.0e-9)      delta = delta_u;
  else if (u <= 0.0+1.0e-9) delta = delta_l;
  else
    {
      if (u <= 0.5)
	{
	  aa = 2.0*u + (1.0-2.0*u)*pow((1+delta_l),(m_distribution + 1.0));
	  delta = pow(aa, (1.0 / (m_distribution + 1.0))) - 1.0;
	}
      else
	{
	  aa = 2.0*(1-u) + 2.0*(u-0.5)*pow((1-delta_u),(m_distribution + 1.0));
	  delta = 1.0 - pow(aa, (1.0 / (m_distribution + 1.0)));
	}
    }
  if(delta < -1.0 || delta > 1.0){
	  throw myException("invalid delta value @ template<typename ED> GAIndividual<ED>::getDelta");
    }
  return (delta);
}

template<typename ED>
 void  GAIndividual<ED>::setDistribution(double dis){
    m_distribution=dis;
 }

 template<typename ED>
void GAIndividual<ED>::setMutationProbability(double p){
    m_mutationProb=p;

}

template<typename ED>
void  GAIndividual<ED>::setXoverProbability(double p){
    m_xoverProb=p;
}

template<typename ED>
void GAIndividual<ED>::realXover(GAIndividual<ED> & mate){
    //Arithmetical crossover
    if(Global::msp_global->mp_uniformAlg->Next()>m_xoverProb) return;
	int numDim=GET_NUM_DIM;
	for(int i=0;i<numDim;i++){	
		if(Global::msp_global->mp_uniformAlg->Next()<=0.5){			
			double r=Global::msp_global->mp_uniformAlg->Next(); //normally r=0.25			
			this->data()[i]= r*this->data()[i]+(1-r)*mate.data()[i];
			mate.data()[i]=r*mate.data()[i]+(1-r)*this->data()[i];						
		}
    }
}

template<typename ED>
void GAIndividual<ED>::singlePointXover(GAIndividual<ED> & mate){
	if(Global::msp_global->mp_uniformAlg->Next()>m_xoverProb) return;
	// single point crossover
	int numDim=GET_NUM_DIM;
	int p=Global::msp_global->getRandInt(0,numDim);

	for (int i = p; i < numDim; i++){
		if(Global::msp_global->mp_uniformAlg->Next()<=0.5){		
			TypeVar x(this->data()[i]);
			this->data()[i]=mate.data()[i];
			mate.data()[i]=x;					  
		}
	}
}

template<typename ED>
void  GAIndividual<ED>::increaseDimension(){
    Individual<ED>::increaseDimension();
    m_backup.increaseDimension();
}

template<typename ED>
void GAIndividual<ED>::decreaseDimension(){
    Individual<ED>::decreaseDimension();
    m_backup.decreaseDimension();

}

#endif // GAINDIVIDUAL_H
