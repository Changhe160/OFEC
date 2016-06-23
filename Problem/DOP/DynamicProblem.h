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
// Created: 11 May 2011
// Last modified:
/// declaration of dynamic problems

#ifndef DYNAMICPROBLEM_H
#define DYNAMICPROBLEM_H

#include "../problem.h"

#define CAST_PROBLEM_DYN dynamic_cast<DynamicProblem*>(Global::msp_global->mp_problem.get())

#define HANDLE_RETURN_FLAG(flag)\
		if(flag==Return_Change_Timelinkage){\
			if(CAST_PROBLEM_DYN) CAST_PROBLEM_DYN->getTriggerTimelinkage()=false;\
		}else if(flag==Return_Change){\
			if(mROOT::getROOT())	mROOT::getROOT()->updateObjValue(Global::msp_global.get());\
		}else if(flag==Return_Terminate) break;

class DynamicProblem : virtual public Problem
{
	struct SChangeType{
		ChangeType type;
		int counter;
		SChangeType & operator=(const SChangeType & rCT){
			if(this==&rCT)  return *this;
			type=rCT.type;
			counter=rCT.counter;
			return *this;
		};
	};
    protected:

        int m_changeFre;
        SChangeType m_changeType;

        int m_period;					// definite period for values repeating
        int m_changeCounter;			// counter of number of changes
        double m_noisySeverity;			// deviation servity from the trajactory of recurrent change
        bool m_flagDimensionChange;			// flag=true, the number of dimensions change, otherwise no change,  default value is false
        bool m_dirDimensionChange;			// direction of change, dir=true means increasing the dimension, otherwise decrease it
        bool m_synchronize;                 // default=true all dimensions change at the same time
		
        int m_dimNumberTemp;                //a temporal variable for dimensional change only

        int m_numPeaks;
        bool m_flagNumPeaksChange;                  // flag of the change of the number of peaks
        bool m_dirNumPeaksChange;                   // true for increasing the number of peaks, otherwise decreasing the number of peaks
        int m_numPeaksTemp;                         // temporal varibal for number of peaks change only


        static const unsigned msc_MaxDimensionNumber=15;
		static const unsigned msc_MinDimensionNumber=2;     //should be greater than 1
	
        static const int msc_MaxNumPeaks=100;
        static const int msc_MinNumPeaks=10;
		
		static thread_local unique_ptr<int> ms_initNumPeaks,ms_initNumDim,ms_numInstance;
        double m_alpha, m_maxAlpha;              // to control step severity
        double m_chaoticConstant;
 
		// features below added on NOV 22 2012
		int m_numPeaksChangeMode;		// for the number of peaks change; 1: periodic with fixed step, 2: periodic with random step, 3: chaotic change
		bool m_noiseFlag;				// objective evaluation with noise in descision space
		bool m_timeLinkageFlag;			// optima move to a random posotion if they are beging tracked 
		
		// severity of noise and time-linkage added in noisy and time-linkage enviroment, June 05 2014
		double m_noiseSeverity_,m_timeLinkageSeverity;
		
		void setDimensionChange(const bool rFlag);
        void setChangeDirction(const bool rFlag);

		bool m_flagTriggerTimeLinkage;
    public:
        static const int msc_NumChangeTypes=12;

        DynamicProblem(const int rId, const int rDimNumber, const int rNumPeaks,const int runId,const unsigned numObj=1 );
        virtual ~DynamicProblem()=0;

        DynamicProblem & operator=(const DynamicProblem & rDP);

        void setChangeFre(const int rChangeFre);
        virtual bool setPeriod(const int rPeriod);
        void setChangeType(const SChangeType &rChangeType);
        void setChangeType(const ChangeType rT);
		void setNumPeaksChange(const bool rPC);
		bool getFlagNumPeaksChange(){
			return m_flagNumPeaksChange;
		}
        void setSynchronize(const bool rFlag);
        void setNoisySeverity(const double rSeverity);
        
        void setAlpha(const double rAlpha){
            m_alpha=rAlpha;
        };
        void setMaxAlpha(const double rMaxAlpha){
            m_maxAlpha=rMaxAlpha;
        };
        void setChoaticConstant(const double rValue){
            m_chaoticConstant=rValue;
        }

        int getChangeFre()const{
            return m_changeFre;
        };
         int getChangeCounter()const {
            return m_changeCounter;
        };
        int getPeriod()const {
            return m_period;
        }
         ChangeType getChangeType() const{
            return m_changeType.type;
        };
         bool getFlagDimensionChange() const{
            return m_flagDimensionChange;
        };
         bool getDirDimensionChange() const{
            return m_dirDimensionChange;
        };
         bool getFlagSynchronizeChange()const{
            return m_synchronize;
        };
    
		void setNumPeakChangeMode(const int mode);
		int getNumPeakChangeMode();
		void setNoiseFlag(const bool flag); 
        int getNumberofPeak()const;
		void setTimeLinkageFlag(const bool flag);
		bool getFlagTimeLinkage(){
			return m_timeLinkageFlag;
		}
        void change();
        double sinValueNoisy(const int x,const double min, const double max, const double amplitude, const double angle,const double noisy_severity=1.);
        double chaoticStep(const double x, const double min, const double max, const double scale=1.0);
		bool predictChange(const int evalsMore);
		void setNoiseSeverity_(double value);
		void setTimeLinkageSeverity(double value);
		bool &getTriggerTimelinkage(){
			return m_flagTriggerTimeLinkage;
		}
		static int getInitialNumPeaks();
    protected:
        virtual void randomChange(){};
        virtual void smallStepChange(){};
        virtual void largeStepChange(){};
        virtual void recurrentChange(){};
        virtual void chaoticChange(){};
        virtual void recurrentNoisyChange(){};

        virtual void changeDimension(){};
        virtual void changeNumPeaks(){};


        void parameterSetting(Problem * rP);
        virtual void  freeMemory(){};
		
};

#endif // DYNAMICPROBLEM_H
