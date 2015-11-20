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
// Last modified:

#include "mMultiModal.h"
#include "../Global/global.h"
#include "mROOT.h"
#include "../Problem/problem.h"
unique_ptr<mMultiModal> mMultiModal::msp_multiModal(nullptr);

mMultiModal::mMultiModal(ParamMap &v):m_avgOptsFound(0), m_avgVisablePeaks(0), m_avgSurvivedPops(0),m_avgInitialRadius(0),
	m_outputProgFlag(false),m_avgPeakQaulity(0),m_avgRadiusQaulity(0),m_avgTrackpercent(0),mpv_infor((int)MAX_NUM_RUN){
    //ctor
	if(Global::msp_global->mp_problem->isProTag(DOP)){
		m_changeFre=(int)(v[param_changeFre]);
	}
  
	 if(!Global::msp_global->mp_problem->isProTag(DOP)){
		m_changeFre=-1;
	}
}

mMultiModal::~mMultiModal()
{
    //dtor
    for(int i=0;MAX_NUM_RUN>i;i++) mpv_infor[i].clear();
	mpv_infor.clear();
}

void mMultiModal::setFileName(stringstream &rName){
    m_fileName.str(rName.str());
}
void mMultiModal::setOutProgFlag(bool rflag){
			m_outputProgFlag=rflag;
}
mMultiModal * mMultiModal::getPopInfor(){
	return mMultiModal::msp_multiModal.get();
}
void mMultiModal::initialize(ParamMap &v){
	mMultiModal::msp_multiModal.reset(new mMultiModal(v));

}
void mMultiModal::deletePerformPopInforDOP(){
	mMultiModal::msp_multiModal.reset();
}
void mMultiModal::input(Global *glob,int fes, int numIndis,int numPops,int numOptsFound,int numVisablePeaks,double initialRadius,\
								double curRadius, double largestRadius, double smallestRadius, double peakQuality,double radiusQaulity, \
								bool isGOptTracked){

	PopInfor infor;
    infor.m_fes=fes;infor.m_numIndis=numIndis;infor.m_numOptsFound=numOptsFound;infor.m_numPops=numPops;infor.m_numVisablePeaks=numVisablePeaks;
    infor.m_avgInitialRadius=initialRadius;infor.m_avgCurRadius=curRadius;infor.m_largestRadius=largestRadius;infor.m_smallestRadius=smallestRadius;
	infor.m_peakQaulity=peakQuality;infor.m_radiusQaulity=radiusQaulity; infor.m_trackPercent=numOptsFound/(double)numVisablePeaks; 
	infor.m_isGOptTracked=isGOptTracked;
	mpv_infor[glob->m_runId].push_back(infor);
}
void mMultiModal::output(){
    /** method to output average results
        example of data format with each input in 3 runs (fes, number of populations), where change fre=10, total fes=50
                    0 10,  0 10,  0 10
                    1 9,   2 9,   1 8
                    3 7,   5 6,   2 7
                    6 2,   7 3,   4 5
                    10 10, 8 1,   6 3
                    12 9,  11 10, 9 1
                    15 6,  13 8,  10 10
                    ...    ...    ...
                     .      .      .
                    49 2   49 1   49 6
    **/
	int maxNumrun=MAX_NUM_RUN;
	int max_fes=mpv_infor[0][mpv_infor[0].size()-1].m_fes;	
	double numSuc=0;
	for(int i=0;i<maxNumrun;i++){
		if(max_fes<mpv_infor[i][mpv_infor[i].size()-1].m_fes)
			max_fes=mpv_infor[i][mpv_infor[i].size()-1].m_fes;
	}
	if(m_changeFre==-1)  //static problem
		m_changeFre=max_fes+1;
	
	vector<unsigned>idx(maxNumrun,0);
  
    int smallestFes;

    stringstream oss;
    ofstream out;

    oss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"PopInfor.txt";
	out.open(oss.str().c_str());
    stringstream ss;
	int lines=0;
	
	while(1){
		int idxSmallest=0;
		while(maxNumrun>idxSmallest&&idx[idxSmallest]>=mpv_infor[idxSmallest].size()) idxSmallest++; // maxNumrun>idxSmallest before the second condition
		if(maxNumrun==idxSmallest){ 
			out<<ss.str();
			break;
		}
		/// find the idex of data with smallest fes
         for(int i=0;maxNumrun>i;i++) {
			 if(idx[i]>=mpv_infor[i].size()) continue;
            if(mpv_infor[i][idx[i]].m_fes<mpv_infor[idxSmallest][idx[idxSmallest]].m_fes) {
                idxSmallest=i;
            }
        }

        smallestFes=mpv_infor[idxSmallest][idx[idxSmallest]].m_fes;
        double numIndis=0, numPops=0, numOptsFound=0, numVisablePeaks=0,initialRadius=0,curRadius=0,largestRadius=0,smallestRadius=0;
        double peakQaulity_=0,radiusQaulity=0,trackPercent=0;
		// sum all the results with smallest idex and the previous results if current index is not the smallest
        for(int i=0;maxNumrun>i;i++){
            if(idx[i]<mpv_infor[i].size()&&mpv_infor[i][idx[i]].m_fes==smallestFes){
                numIndis+=mpv_infor[i][idx[i]].m_numIndis;
                numPops+=mpv_infor[i][idx[i]].m_numPops;
                numOptsFound+=mpv_infor[i][idx[i]].m_numOptsFound;
                numVisablePeaks+=mpv_infor[i][idx[i]].m_numVisablePeaks;
                initialRadius+=mpv_infor[i][idx[i]].m_avgInitialRadius;
                curRadius+=mpv_infor[i][idx[i]].m_avgCurRadius;
                largestRadius+=mpv_infor[i][idx[i]].m_largestRadius;
                smallestRadius+=mpv_infor[i][idx[i]].m_smallestRadius;
				peakQaulity_+=mpv_infor[i][idx[i]].m_peakQaulity;
				radiusQaulity+=mpv_infor[i][idx[i]].m_radiusQaulity;
				trackPercent+=mpv_infor[i][idx[i]].m_trackPercent;
                idx[i]++;

            }else{
				int indexPre=0;
				if(idx[i]>0) indexPre=idx[i]-1;
				
				numIndis+=mpv_infor[i][indexPre].m_numIndis;
				numPops+=mpv_infor[i][indexPre].m_numPops;
				numOptsFound+=mpv_infor[i][indexPre].m_numOptsFound;
				numVisablePeaks+=mpv_infor[i][indexPre].m_numVisablePeaks;
				initialRadius+=mpv_infor[i][indexPre].m_avgInitialRadius;
                curRadius+=mpv_infor[i][indexPre].m_avgCurRadius;
                largestRadius+=mpv_infor[i][indexPre].m_largestRadius;
                smallestRadius+=mpv_infor[i][indexPre].m_smallestRadius;
				peakQaulity_+=mpv_infor[i][indexPre].m_peakQaulity;
				radiusQaulity+=mpv_infor[i][indexPre].m_radiusQaulity;
				trackPercent+=mpv_infor[i][indexPre].m_trackPercent;
            }
        }
		bool flag=true;
		for(int i=0;maxNumrun>i;i++){
			if(idx[i]<1){ flag=false; break;}

			if(idx[i]==mpv_infor[i].size()) continue; //make sure the last data can be used

			if(mpv_infor[i][idx[i]].m_fes/m_changeFre==mpv_infor[i][idx[i]-1].m_fes/m_changeFre)
			{
				flag=false; break;
			}
		}
		if(flag){
			for(int i=0;maxNumrun>i;i++) if(mpv_infor[i][idx[i]-1].m_isGOptTracked) numSuc+=1;
            m_avgOptsFound+=numOptsFound;
            m_avgVisablePeaks+=numVisablePeaks;
            m_avgSurvivedPops+=numPops;
            m_avgInitialRadius+=initialRadius;
			m_avgPeakQaulity+=peakQaulity_;
			m_avgRadiusQaulity+=radiusQaulity;
			m_avgTrackpercent+=trackPercent;
		}

        numIndis/=maxNumrun;
        numPops/=maxNumrun;
        numOptsFound/=maxNumrun;
        numVisablePeaks/=maxNumrun;
        initialRadius/=maxNumrun;
        curRadius/=maxNumrun;
        smallestRadius/=maxNumrun;
        largestRadius/=maxNumrun;
		trackPercent/=maxNumrun;
        if(m_outputProgFlag){
			//ss<<smallestFes<<" "<<numIndis<<" "<<numPops<<" "<<numOptsFound<<" "<<trackPercent<<" "<<numVisablePeaks<<" "<<initialRadius<<" "<<curRadius<<" "<<smallestRadius<<" "<<largestRadius<<endl;
			ss<<smallestFes<<" "<<numIndis<<" "<<numPops<<" "<<numOptsFound<<" "<<numVisablePeaks<<endl;
			
			lines++;
			if(lines%2000==0){
				out<<ss.str();
				ss.str("");
			}

		}
	}
	out.close();

	// calculate the average tracking percent and average tracking rate for the GOpt
	//double **avg_per,**avg_gRatio;
	

	unsigned numChanges=max_fes/m_changeFre+1;
	vector< vector<double> > avg_per(numChanges,vector<double>(maxNumrun+2));
	vector< vector<double> > avg_gRatio(numChanges,vector<double>(maxNumrun+2));
	
	
	for(int j=0;maxNumrun>j;j++){
		unsigned k=1;	
		for(unsigned i=0;i<numChanges;i++){
			if(mpv_infor[j].size()<numChanges) throw myException("missing data@mMultiModal::output()");
			while(k<mpv_infor[j].size()&&mpv_infor[j][k].m_fes/m_changeFre==mpv_infor[j][k-1].m_fes/m_changeFre) k++;
	
			avg_per[i][j]=mpv_infor[j][k-1].m_trackPercent;
			if(mpv_infor[j][k-1].m_isGOptTracked)	avg_gRatio[i][j]=1.;
			else avg_gRatio[i][j]=0;
			k++;
		}
		
	}

	for(unsigned i=0;i<numChanges;i++){
		avg_gRatio[i][maxNumrun]=0;
		for(int j=0;maxNumrun>j;j++){
			avg_gRatio[i][maxNumrun]+=avg_gRatio[i][j];
		}
		avg_gRatio[i][maxNumrun]/=maxNumrun;
		avg_gRatio[i][maxNumrun+1]=0;
		for(int j=0;maxNumrun>j;j++){
			avg_gRatio[i][maxNumrun+1]+=(avg_gRatio[i][j]-avg_gRatio[i][maxNumrun])*(avg_gRatio[i][j]-avg_gRatio[i][maxNumrun]);
		}
		avg_gRatio[i][maxNumrun+1]=sqrt(avg_gRatio[i][maxNumrun+1]/maxNumrun);
	}

	for(unsigned i=0;i<numChanges;i++){
		avg_per[i][maxNumrun]=0;
		for(int j=0;maxNumrun>j;j++){
			avg_per[i][maxNumrun]+=avg_per[i][j];
		}
		avg_per[i][maxNumrun]/=maxNumrun;
		avg_per[i][maxNumrun+1]=0;
		for(int j=0;maxNumrun>j;j++){
			avg_per[i][maxNumrun+1]+=(avg_per[i][j]-avg_per[i][maxNumrun])*(avg_per[i][j]-avg_per[i][maxNumrun]);
		}
		avg_per[i][maxNumrun+1]=sqrt(avg_per[i][maxNumrun+1]/maxNumrun);
	}
	double var1=0,var2=0,avgPercent=0,avgGoptRatio=0;
	for(unsigned i=0;i<numChanges;i++){
		avgPercent+=avg_per[i][maxNumrun];
		var1+=avg_per[i][maxNumrun+1];
		avgGoptRatio+=avg_gRatio[i][maxNumrun];
		var2+=avg_gRatio[i][maxNumrun+1];
	}
	avgPercent/=numChanges;var1/=numChanges;
	avgGoptRatio/=numChanges;var2/=numChanges;

	oss.str("");

	oss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Sta.txt";
	out.open(oss.str().c_str(),ios_base::app);
    out<<"peaksFound: "<<m_avgOptsFound/(maxNumrun*(numChanges))<<endl;
	out<<"success rate: "<<numSuc/(maxNumrun*(numChanges))<<endl;
    //out<<"trackingPercent: "<<m_avgTrackpercent/(maxNumrun*(numChanges))<<endl;
	out<<"peakRatio&Var: "<<avgPercent<<" "<<var1<<endl;
	//out<<"TrackingRatioGopt&Var: "<<avgGoptRatio<<" "<<var2<<endl;
	//out<<"SurvivedPops: "<<m_avgSurvivedPops/(maxNumrun*(numChanges))<<endl;
    out<<"VisablePeaks: "<<m_avgVisablePeaks/(maxNumrun*(numChanges))<<endl;
   // out<<"initialRadius: "<<m_avgInitialRadius/(maxNumrun*(numChanges))<<endl;
	out<<"PeaksFoundQaulity: "<<m_avgPeakQaulity/(maxNumrun*(numChanges))<<endl;
	//out<<"RadiusQaulity: "<<m_avgRadiusQaulity/(maxNumrun*(numChanges))<<endl;
	//ROOT::getROOT()->calculateMean();
	//out<<"ROOT&var: "<<ROOT::getROOT()->getMeanRoot()<<" "<<ROOT::getROOT()->getROOTStd()<<endl;
	out.close();
}
void mMultiModal::setFileName(ParamMap &v){
	m_fileName.str("");	
	for(auto &i:v){
		for(auto &j:Global::msm_param){
			if(i.first==param_gOptFlag||i.first==param_algId||i.first==param_proId||i.first==param_flagNoise||\
				i.first==param_flagNumPeakChange||i.first==param_flagTimeLinkage||i.first==param_numRun||\
				i.first==param_numTask||i.first==param_minNumPopSize||i.first==param_hibernatingRadius||\
				i.first==param_solutionValidationMode||i.first==param_evalCountFlag||\
				i.first==param_workingDir||i.first==param_sampleFre||i.first==param_maxEvals||i.first==param_flagNumPeakChange||\
				i.first==param_peakNumChangeMode) continue;
			if(i.first==j.second){			
				m_fileName<<j.first.substr(6)<<i.second<<"_";			
				break;
			}
		}		
	}
}

void mMultiModal::deleteMultiModal(){
	msp_multiModal.reset();
}