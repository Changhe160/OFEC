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
// Created: 11 May 2011
// Last modified:
#ifndef HYBRIDCOMP_H
#define HYBRIDCOMP_H

#include "BenchmarkFunction.h"

class HybridComp : public BenchmarkFunction
{
    public:
        HybridComp(ParamMap &v);
        HybridComp(const int rId, const int rDimNumber, string& rName);
        virtual ~HybridComp();
        void setFunction(unsigned * rId, string rFucs[]);
        int getNumFuncs();
   
    protected:

        virtual void allocateMemory(const int rNumDim, const int rNumFuc);
        virtual void freeMemory();
        virtual void initialize();
        virtual bool loadTranslation();
        virtual bool loadRotation();

        void setUpFCom();
        void SetUpFCom_CEC05();
        void SetUpFRH_Com_CEC05();
        void SetUpFRH001_Com_CEC05();
        void SetUpFRH002_Com_CEC05();
		void evaluate__(double const  *x,vector<double>& obj);

    private:
        static const  int  m_numFuncs=10;				// number of basic functions, for hybrid functions
        BenchmarkFunction **mpp_f; // the functions

        double *mp_convergeSeverity;        // severity of converge range for each function
        double *mp_stretchSeverity;          // severity of stretching original function, greater than 1 for stretch
        double *mp_weight;                  // weight value of each basic function                            
        double *mp_height;
        double *mp_fmax;

        double m_heightNormalizeSeverity;   // constant number for noralizing all basic function with similar height
};

#endif // HYBRIDCOMP_H
