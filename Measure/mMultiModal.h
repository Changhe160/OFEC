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
// Created: 28 September 2011
// Last modified:
/*
    This class is used to record relevant inter-information for algorithms using multiple population methods in DOPs during run time
*/

#ifndef POPULATIONINFORDOP_H
#define POPULATIONINFORDOP_H
class Global;
#include "../Utility/include.h"
#include "../Utility/TypeVar/typeVar.h"

class mMultiModal
{
    struct PopInfor{
        int m_fes, m_numIndis,m_numPops,m_numOptsFound,m_numVisablePeaks;
        double m_avgInitialRadius,m_avgCurRadius,m_largestRadius,m_smallestRadius;
		double m_peakQaulity,m_radiusQaulity;
		
		//12/07/2013
		double m_trackPercent; // the percentage of tracked peaks over all peaks
		bool m_isGOptTracked;
    };
    private:
        mMultiModal(ParamMap &v);
        vector<vector<PopInfor>> mpv_infor;
        stringstream m_fileName;
		bool m_outputProgFlag;

		double m_avgOptsFound, m_avgVisablePeaks, m_avgSurvivedPops,m_avgInitialRadius;
		double m_avgPeakQaulity,m_avgRadiusQaulity;
		double m_avgTrackpercent;
		unsigned m_changeFre;
		static unique_ptr<mMultiModal>  msp_multiModal;
    public:

		void input(Global*, int fes, int numIndis,int numPops,int numOptsFound,int numVisablePeaks,double initialRadius=0,double curRadius=0,double largestRadius=0,\
			double smallestRadius=0, double peakQaulity=0,double radiusQaulity=0, bool isGOptTracked=false);
        void output();
        void setFileName(stringstream &rName);
		void setOutProgFlag(bool rflag);
        virtual ~mMultiModal();
       
        static mMultiModal * getPopInfor();
        static void initialize(ParamMap &v);
        static void deletePerformPopInforDOP();
		void setFileName(ParamMap &v);
		static void deleteMultiModal();
};

#endif // POPULATIONINFORDOP_H
