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

#ifndef DYNAMICCONTINUOUS_H
#define DYNAMICCONTINUOUS_H

#include "../ContinuousProblem.h"
#include "DynamicProblem.h"

#define CAST_PROBLEM_DYN_CONT dynamic_cast<DynamicContinuous*>(Global::msp_global->mp_problem.get())

class DynamicContinuous :  public DynamicProblem, public ContinuousProblem
{
    protected:
 //       int m_numPeaks;						        // number of peaks in Rotation_DBG , number of function in Composition_DBG
        double **mpp_peak;					    	// positions of local or global optima(local optima in Rotation_DBG,
        double **mpp_prePeak;							// global optima of basic function in Composition_DBG)
        double **mpp_initialPeak;				        // save the initial positions
        double *mp_height;							// peak height in Rotation_DBG, height of global optima in Composition_DBG
        double *mp_width;                           // weight value of each basic function in Composition_DBG,  peak width in Rotation_DBG

        ///TODO preHeight and preWidth not considered in current version
        double *mp_preHeight;
        double *mp_preWidth;

        double m_minHeight,m_maxHeight;		// minimum\maximum height of all peaks(local optima) in Rotation_DBG(Composition_DBG)
        double *mp_heightSeverity;

        double m_minWidth, m_maxWidth;
        double *mp_widthSeverity;

        double *mp_fit;						    	// objective value of each basic funciton in Composition_DBG, peak height in Rotation_DBG
        double m_globalOptima;				    	// global optima value
        bool *mp_globalOptimaIdx;                      // the index of the global optimal peak

        int m_currentPeak;                         // the peak where the best individual is located
        int m_maxPeaksNumber;                      // the number of heigthest peaks
        double m_currentBest;                      // the objective value of current best individual

        bool *mp_whetherChange;                      // whether peaks change or not
        int m_numChangePeaks;                       // the number of peaks that change
        double m_changePeakRatio;                    // the ratio of changing peaks

        int m_numVisablePeaks;                      // number of visable peaks, a peak is visable only if no peak is on top of it
		int *mp_isTracked;							// accumulated number of peak[i] being tracked      
		int *mp_heightOrder;
		int m_peaksFound;
		bool *mp_found;	
		double *mp_timeLinkage;						//a random vector to change peaks position which are being tracked
		int *mp_amendedHeightOrder;
	   //added 16/05/2013
		double *mp_associateRadius; /*actual radius in the fitness landscape*/
		double m_peakQaulity;	/*to evaluate qaulity of peaks trcked  in terms of peaks heights*/
		//added 04/07/2014
		bool isGOpt(int idx);
public:
        DynamicContinuous(const int rId, const int rDimNumber,  const int rNumPeaks,const int runId,const unsigned numObj=1);
        virtual ~DynamicContinuous()=0;

        double getGlobalMax()const;
        void printFun( std::ofstream & out);
        const double * getPeak(const int p) const;
        const double * const* const getAllPeaks()const;
        double getPeakHeight(const int p)const;
        double getPrePeakHeight(const int p)const;
        double getPrePeakWidth(const int p)const;
        const double *const getHeight() const;
        const double *const getPrePeak(const int p)const;
        const bool *const getGlobalOptimaIdx()const;
        int getNumberofGlobalOptPeak()const;

        void setNumberofChanges(const int n);
        void setNumberofChanges(const double rRatio);

        void setHeightSeverity(const double rS);
        void setWidthSeverity(const double rS);
        DynamicContinuous &operator=(const DynamicContinuous &rDCP);

        void setHeight(const double *h);
        void setPosition(const double * const * const p);//const double **p
        virtual void setWidth(const double w);

        /// for debug mode
        void printPeak( const int rIdx);
        void printPeaks(ofstream & out);
        int getNumofVisablePeaks();
        bool isVisable(const int rIdx);
		int getTrackNumber(int idex);
		bool isTracked(vector<double> &gen,vector<double> &obj); // is any peak tracked for the first time
		bool isTracked(double *gen,vector<double> &obj);
		int getPeaksFound();
		double getAssociateRadius( int idx);
		double getPeaksTracedQaulity();
		//15-07-2013
		bool isGOptTracked();
		const double * getNearestPeak(const vector<double>& );
    protected:
        virtual void randomChange(){};
        virtual void smallStepChange(){};
        virtual void largeStepChange(){};
        virtual void recurrentChange(){};
        virtual void chaoticChange(){};
        virtual void recurrentNoisyChange(){};

        void parameterSetting(Problem * rP);
        virtual void  freeMemory();
        virtual void allocateMemory(const int rDimNum, const int rPeaks);

        virtual void changeDimension(){};
        virtual void changeNumPeaks(){};

        void calculateGlobalOptima();
        void updateNumberofChanges();
        void computeNumVisablePeaks();
		void addNoise(double *x);
		void updateTimeLinkage();
		void movePeak(const int idx);
		virtual void calculateAssociateRadius(){};
		virtual void updatePeakQaulity(){};

};

#endif // DYNAMICCONTINUOUS_H
