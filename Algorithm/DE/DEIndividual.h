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
// Created: 21 September 2011
// Last modified: 12 Dec. 2014
#ifndef DEINDIVIDUAL_H
#define DEINDIVIDUAL_H
#include "../Individual.h"

template<class,  class> class DEPopulation;
class DEIndividual: public Individual<CodeVReal>
{
	template<class ,class> friend class DEPopulation;
protected:
        Solution<CodeVReal> m_pv,m_pu;    // donor vector and trial vector, respectively.
    public:
		 DEIndividual();
        virtual ~DEIndividual();
		DEIndividual( const DEIndividual &p); 
		DEIndividual( const Solution<CodeVReal> &p); 
		DEIndividual & operator=(const DEIndividual &other);
        virtual ReturnFlag initialize(bool mode=true);
        virtual ReturnFlag initialize(int idex, int id,bool mode=true);
        virtual void initialize(const Solution<CodeVReal> &p, int idex, int id);
		virtual ReturnFlag initialize(int rIdx,int rID, int rPopsize,bool mode=true );
		virtual ReturnFlag initialize(const Solution<CodeVReal> &p,double radius,  int ridex, int rid,bool mode=true);
		void increaseDimension();
		void decreaseDimension();
        virtual void mutate(double F, Solution<CodeVReal> *r1,  Solution<CodeVReal> *r2,  Solution<CodeVReal> *r3, Solution<CodeVReal> *r4=0,  Solution<CodeVReal> *r5=0);
        virtual void recombine(double CR);
        virtual ReturnFlag select();
		void printToScreen();
};

#endif // DEINDIVIDUAL_H
