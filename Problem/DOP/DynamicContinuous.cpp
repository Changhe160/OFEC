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

#include "DynamicContinuous.h"
#include "../../Global/global.h"
#include "../../Utility/myVector.h"
#ifdef DEMON_OFEC
#include <windows.h>
#include <process.h>   
extern void calculateSamplePoints();
extern HANDLE g_hthDraw;
#endif

DynamicContinuous::DynamicContinuous(const int rId, const int rDimNumber,  const int rNumPeaks,const int runId,const unsigned numObj):\
	DynamicProblem(rId,rDimNumber, rNumPeaks,numObj),ContinuousProblem(rId,rDimNumber,string(),numObj),\
	m_numChangePeaks(m_numPeaks),m_changePeakRatio(1.0){
	setAccuracy(0.01);
	allocateMemory(m_numDim,m_numPeaks);
	for (int i=0; i< m_numPeaks; i++) mp_whetherChange[i]=true;
    m_proPar<<"Changing peaks ratio:"<<m_changePeakRatio<<"; ";
}

DynamicContinuous::~DynamicContinuous(){
    freeMemory();
}

void DynamicContinuous::allocateMemory(const int rDimNum, const int rPeaks){

	mpp_peak = new double*[rPeaks];
	mpp_prePeak= new double*[rPeaks];
    mpp_initialPeak=new double*[rPeaks];
    mp_width=new double[rPeaks];
    mp_height=new double[rPeaks];
    mp_preHeight=new double[rPeaks];
    mp_preWidth=new double[rPeaks];
    mp_fit=new double[rPeaks];
	mp_widthSeverity=new double [rPeaks];
	mp_heightSeverity=new double[rPeaks];

	for (int i=0; i< rPeaks; i++){
        mpp_peak[i]= new double[rDimNum];
        mpp_prePeak[i]= new double[rDimNum];
        mpp_initialPeak[i]= new double[rDimNum];
	}
	mp_whetherChange=new bool[rPeaks];
	mp_globalOptimaIdx=new bool[rPeaks];

	mp_isTracked= new int [rPeaks];
	mp_heightOrder=new int [rPeaks];
	mp_found=new bool[rPeaks];
	
	mp_timeLinkage=new double[rPeaks];
	mp_amendedHeightOrder=new int[rPeaks];
	mp_associateRadius = new double[rPeaks];
	for(int i=0;i<rPeaks;++i){
		mp_found[i]=false;
		mp_globalOptimaIdx[i]=false;
		mp_whetherChange[i]=true;
		mp_associateRadius[i]=0;
		mp_heightOrder[i]=-1;
	} 
}
void  DynamicContinuous::freeMemory(){
	for (int i=0; i< m_numPeaks; i++){
		delete [] mpp_peak[i];
		delete [] mpp_prePeak[i];
		delete[] mpp_initialPeak[i];

	}
	delete [] mpp_peak;
	delete [] mpp_prePeak;
	delete[] mpp_initialPeak;
	delete[] mp_height;
	delete[] mp_width;
	delete[] mp_preHeight;
	delete[] mp_preWidth;
	delete[] mp_fit;
	delete [] mp_whetherChange;
	delete [] mp_globalOptimaIdx;

	delete [] mp_isTracked;
	delete [] mp_heightOrder;
	delete [] mp_found;
	delete []mp_timeLinkage;
	delete [] mp_amendedHeightOrder;
	mp_isTracked=0;
	mp_heightOrder=0;
	mp_found=0;

    mpp_peak=0;
    mpp_prePeak=0;		mpp_initialPeak=0;
    mp_height=0;		mp_width=0;
    mp_preHeight=0;		mp_preWidth=0;
    mp_fit=0;			mp_whetherChange=0;
	mp_timeLinkage=0;
	mp_amendedHeightOrder=0;
	delete [] mp_associateRadius;
	mp_associateRadius=0;
	delete [] mp_widthSeverity;
	delete [] mp_heightSeverity;
	mp_widthSeverity=0;
	mp_heightSeverity=0;
}

double DynamicContinuous::getGlobalMax()const{
	return m_globalOptima;
}

 void DynamicContinuous::printFun(ofstream & out){

	 for(int i=0;i<m_numPeaks;i++){

		 for(int j=0;j<m_numDim;j++)
			 out<<mpp_peak[i][j]<<" ";
		out<<endl;
	 }
 }

const double * DynamicContinuous::getPeak(const int p)const {
	 if(p<0||p>=m_numPeaks) {
		 throw myException("Please give right value of peak index [0,] @DynamicContinuous::getPeak");
	 }
	return mpp_peak[p];
}

const double * const DynamicContinuous::getPrePeak(const int p) const{
	 if(p<0||p>=m_numPeaks) {
		 throw myException("Please give right value of peak index [0,] @DynamicContinuous::getPrePeak ");
	 }
	return  mpp_prePeak[p];
}
const double * const * const DynamicContinuous::getAllPeaks()const {
	 return  mpp_peak;
 }

 double DynamicContinuous::getPeakHeight(const int p)const{

	if(p<0||p>=m_numPeaks) {
		 throw myException("Please give right value of peak index [0,] @DynamicContinuous::getPeakHeight");
	 }
    return mp_height[p];
}


 double DynamicContinuous::getPrePeakHeight(const int p)const {

    if(p<0||p>=m_numPeaks) {
		 throw myException("Please give right value of peak index [0,] @DynamicContinuous::getPrePeakHeight");
	 }
    return mp_preHeight[p];

}
double DynamicContinuous::getPrePeakWidth(const int p)const{

     if(p<0||p>=m_numPeaks) {
		 throw myException("Please give right value of peak index [0,] DynamicContinuous::getPrePeakWidth");
	 }
    return mp_preWidth[p];
}
const bool * const DynamicContinuous::getGlobalOptimaIdx() const{

	 return mp_globalOptimaIdx;
 }
int DynamicContinuous::getNumberofGlobalOptPeak()const{
    return m_maxPeaksNumber;
}

 void DynamicContinuous::setNumberofChanges(const int n){
	 if(n<1||n>m_numPeaks) {
        throw myException("the number of changing peaks is invalid@DynamicContinuous::setNumberofChanges");
	 }

	m_numChangePeaks=n;
	m_changePeakRatio=n/(double)m_numPeaks;
	updateNumberofChanges();

    size_t start, end;
    start=m_proPar.str().find("Changing peaks ratio:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Changing peaks ratio:"<<m_changePeakRatio<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);

 }
void DynamicContinuous::setNumberofChanges(const double rRatio){
	 if(rRatio<0||rRatio>1) {
        throw myException("the ratio of changing peaks is invalid@DynamicContinuous::setNumberofChanges");
	 }

	m_changePeakRatio=rRatio;
	m_numChangePeaks=(int)(m_numPeaks*m_changePeakRatio)>1?(int)(m_numPeaks*m_changePeakRatio):1;
    updateNumberofChanges();

    size_t start, end;
    start=m_proPar.str().find("Changing peaks ratio:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Changing peaks ratio:"<<m_changePeakRatio<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);

 }

 void DynamicContinuous::updateNumberofChanges(){
	 if(m_numChangePeaks==m_numPeaks){
	 for(int i=0;i<m_numPeaks;i++) mp_whetherChange[i]=true;
	 return;
	 }
	int *a=new int[m_numPeaks];

	Global::msp_global->initializeRandomArray(a,m_numPeaks,Program_Problem);
	// make sure the global optimum changes always
	int gopt=0;
	for(int i=0;i<m_numPeaks;i++){
	     if( mp_globalOptimaIdx[i]){
			 gopt=i;
			 break;
		 }
	}
	int gidx;
	for(int i=0;i<m_numPeaks;i++){
		if( a[i]==gopt){
			 gidx=i;
			 break;
		 }
	}
	int t=a[0];
	a[0]=a[gidx];
	a[gidx]=t;

	for(int i=0;i<m_numPeaks;i++) mp_whetherChange[i]=false;
	for(int i=0;i<m_numChangePeaks;i++) mp_whetherChange[a[i]]=true;
	delete []a;

 }
void DynamicContinuous::setHeightSeverity(const double rS){
    for(int i=0;i<m_numPeaks;i++) 	mp_heightSeverity[i]=rS;
}
void DynamicContinuous::setWidthSeverity(const double rS){
    for(int i=0;i<m_numPeaks;i++) 	mp_widthSeverity[i]=rS;
}
DynamicContinuous & DynamicContinuous::operator=(const DynamicContinuous &rDCP){
    if(this==&rDCP) return *this;

    if(m_numDim!=rDCP.m_numDim||m_numPeaks!=rDCP.m_numPeaks) return *this;

    DynamicProblem::operator=(rDCP);
	ContinuousProblem::operator=(rDCP);

    for(int i=0;i<m_numPeaks;i++){
        copy(rDCP.mpp_peak[i],rDCP.mpp_peak[i]+m_numDim,mpp_peak[i]);
        copy(rDCP.mpp_prePeak[i],rDCP.mpp_prePeak[i]+m_numDim,mpp_prePeak[i]);
        copy(rDCP.mpp_initialPeak[i],rDCP.mpp_initialPeak[i]+m_numDim,mpp_initialPeak[i]);
    }
    copy(rDCP.mp_height,rDCP.mp_height+m_numPeaks,mp_height);
    copy(rDCP.mp_width,rDCP.mp_width+m_numPeaks,mp_width);
    copy(rDCP.mp_preHeight,rDCP.mp_preHeight+m_numPeaks,mp_preHeight);
    copy(rDCP.mp_preWidth,rDCP.mp_preWidth+m_numPeaks,mp_preWidth);
    copy(rDCP.mp_fit,rDCP.mp_fit+m_numPeaks,mp_fit);
    copy(rDCP.mp_whetherChange,rDCP.mp_whetherChange+m_numPeaks,mp_whetherChange);

    m_minHeight=rDCP.m_minHeight;
    m_maxHeight=rDCP.m_maxHeight;
   
    m_minWidth=rDCP.m_minWidth;
    m_maxWidth=rDCP.m_maxWidth;
   
	copy(rDCP.mp_heightSeverity,rDCP.mp_heightSeverity+m_numPeaks,mp_heightSeverity);
	copy(rDCP.mp_widthSeverity,rDCP.mp_widthSeverity+m_numPeaks,mp_widthSeverity);

    m_globalOptima=rDCP.m_globalOptima;
    copy(rDCP.mp_globalOptimaIdx,rDCP.mp_globalOptimaIdx+m_numPeaks,mp_globalOptimaIdx);

    m_currentBest=rDCP.m_currentBest;
    m_currentPeak=rDCP.m_currentPeak;
    m_maxPeaksNumber=rDCP.m_maxPeaksNumber;

    m_numChangePeaks =rDCP.m_numChangePeaks;
    m_changePeakRatio=rDCP.m_changePeakRatio;

	m_numVisablePeaks=rDCP.m_numVisablePeaks;

	copy(rDCP.mp_isTracked,rDCP.mp_isTracked+m_numPeaks,mp_isTracked);
	copy(rDCP.mp_heightOrder,rDCP.mp_heightOrder+m_numPeaks,mp_heightOrder);
	copy(rDCP.mp_found,rDCP.mp_found+m_numPeaks,mp_found);
	m_peaksFound=rDCP.m_peaksFound;
	copy(rDCP.mp_timeLinkage,rDCP.mp_timeLinkage+m_numPeaks,mp_timeLinkage);
	
	m_searchRange=rDCP.m_searchRange;

	copy(rDCP.mp_amendedHeightOrder,rDCP.mp_amendedHeightOrder+m_numPeaks,mp_amendedHeightOrder);
	copy(rDCP.mp_associateRadius,rDCP.mp_associateRadius+m_numPeaks,mp_associateRadius);
	m_peakQaulity=rDCP.m_peakQaulity;
	m_globalOpt=rDCP.m_globalOpt;
	
    return *this;
}
void DynamicContinuous::parameterSetting(Problem * rP){
    DynamicProblem::parameterSetting(rP);
	ContinuousProblem::parameterSetting(rP);

	DynamicContinuous *dcp=dynamic_cast<DynamicContinuous *>(rP);

	int dim=m_dimNumberTemp<rP->getNumDim()?m_dimNumberTemp:rP->getNumDim();
    int peaks=m_numPeaks<dcp->getNumberofPeak()?m_numPeaks:dcp->getNumberofPeak();

    for(int i=0;i<peaks;i++){
   		copy(dcp->mpp_peak[i],dcp->mpp_peak[i]+m_numDim,mpp_peak[i]);
        copy(dcp->mpp_prePeak[i],dcp->mpp_prePeak[i]+m_numDim,mpp_prePeak[i]);
        copy(dcp->mpp_initialPeak[i],dcp->mpp_initialPeak[i]+m_numDim,mpp_initialPeak[i]);
    }
    copy(dcp->mp_height,dcp->mp_height+peaks,mp_height);
    copy(dcp->mp_width,dcp->mp_width+peaks,mp_width);
    copy(dcp->mp_preHeight,dcp->mp_preHeight+peaks,mp_preHeight);
    copy(dcp->mp_preWidth,dcp->mp_preWidth+peaks,mp_preWidth);
    copy(dcp->mp_fit,dcp->mp_fit+peaks,mp_fit);
    copy(dcp->mp_whetherChange,dcp->mp_whetherChange+peaks,mp_whetherChange);

    m_minHeight=dcp->m_minHeight;
    m_maxHeight=dcp->m_maxHeight;

    m_minWidth=dcp->m_minWidth;
    m_maxWidth=dcp->m_maxWidth;
    
	copy(dcp->mp_heightSeverity,dcp->mp_heightSeverity+peaks,mp_heightSeverity);
	copy(dcp->mp_widthSeverity,dcp->mp_widthSeverity+peaks,mp_widthSeverity);

    m_globalOptima=dcp->m_globalOptima;
	copy(dcp->mp_globalOptimaIdx,dcp->mp_globalOptimaIdx+peaks,mp_globalOptimaIdx);
    
    m_currentBest=dcp->m_currentBest;
    m_currentPeak=dcp->m_currentPeak;
    m_maxPeaksNumber=dcp->m_maxPeaksNumber;

    m_changePeakRatio=dcp->m_changePeakRatio;
    m_numChangePeaks =(int)(m_changePeakRatio*peaks)>1?(int)(m_changePeakRatio*peaks):1;//dcp->m_numChangePeaks;

	copy(dcp->mp_isTracked,dcp->mp_isTracked+peaks,mp_isTracked);
	copy(dcp->mp_heightOrder,dcp->mp_heightOrder+peaks,mp_heightOrder);
	copy(dcp->mp_found,dcp->mp_found+peaks,mp_found);
	m_peaksFound=dcp->m_peaksFound;
	copy(dcp->mp_timeLinkage,dcp->mp_timeLinkage+peaks,mp_timeLinkage);
	copy(dcp->mp_amendedHeightOrder,dcp->mp_amendedHeightOrder+peaks,mp_amendedHeightOrder);
	copy(dcp->mp_associateRadius,dcp->mp_associateRadius+peaks,mp_associateRadius);

	for(int j=0;j<dim;j++){
        m_searchRange[j]=dcp->m_searchRange[j];
    }

    m_peakQaulity=dcp->m_peakQaulity;

}

void DynamicContinuous::calculateGlobalOptima(){

	if(m_OptMode[0]==MAX_OPT) m_globalOptima=*max_element(mp_height,mp_height+m_numPeaks);
	else m_globalOptima=*min_element(mp_height,mp_height+m_numPeaks);
	
	m_globalOpt.clear();
    m_maxPeaksNumber=0;
	double mindis=LONG_MAX;
	for(int i=0;i<m_numPeaks;i++){
	    mp_globalOptimaIdx[i]=false;	
		if(mp_height[i]==m_globalOptima){
			for(int j=0;j<m_numPeaks;++j){
				if(j==i) continue;
				MyVector s1(m_numDim,mpp_peak[i]),s2(m_numDim,mpp_peak[j]);

				double dis=s1.getDis(s2);
				if(mindis>dis){
					mindis=dis;
				}
			}
		    m_maxPeaksNumber++;
		    mp_globalOptimaIdx[i]=true;
			Solution<CodeVReal> s(m_numDim,m_numObj);
			copy(mpp_peak[i],mpp_peak[i]+m_numDim,s.data().m_x.begin());
			s.data().m_obj[0]=mp_height[i];
			m_globalOpt.appendOptima(s,false);
		}
	}

	if(IS_PROBLEM_NAME(m_id,"DYN_CONT_RotationDBG")||IS_PROBLEM_NAME(m_id,"DYN_CONT_MovingPeak")){
		if(mindis/2<m_disAccuracy)		setDisAccuracy(mindis/2);
	}
    computeNumVisablePeaks();

	m_peaksFound=0;
	if(m_timeLinkageFlag) updateTimeLinkage();
	for (int i=0; i< m_numPeaks; i++){
		mp_heightOrder[i]=i;
		mp_found[i]=false;
	}
	vector<int> idx(m_numPeaks);
	gQuickSort(mp_height,m_numPeaks,idx);
	copy(idx.begin(),idx.end(),mp_heightOrder);
	gAmendSortedOrder<double*>(mp_height,mp_heightOrder,mp_amendedHeightOrder,m_numPeaks);
	//calculateAssociateRadius();
	m_peakQaulity=0;
}

void DynamicContinuous::setHeight(const double *h){
    copy(h,h+m_numPeaks,mp_height);
}


//const double **p

void DynamicContinuous::setPosition(const double * const * const p){
    for(int i=0;i<m_numPeaks;i++){
        copy(p[i],p[i]+m_numDim,mpp_peak[i]);
        copy(p[i],p[i]+m_numDim,mpp_initialPeak[i]);
    }
}
const double *const DynamicContinuous::getHeight() const{
    return mp_height;
}
void DynamicContinuous::setWidth(const double w){
    for(int i=0;i<m_numPeaks;i++)
        mp_width[i]=w;
}
void DynamicContinuous::printPeak( const int rIdx){

    cout<<"the "<<rIdx<<"th peak, height: "<<mp_height[rIdx]<<" ass radius: "<<mp_associateRadius[rIdx]<<" position:"<<endl;

    for(int i=0;i<m_numDim;i++){
        cout<<mpp_peak[rIdx][i]<<" ";
    }
    cout<<endl;
}
void DynamicContinuous::printPeaks(ofstream & out){
    for(int j=0;j<m_numPeaks;j++){
         for(int i=0;i<m_numDim;i++){
            out<<mpp_peak[j][i]<<" ";
        }
        out<<mp_height[j]<<endl;
    }

	/* for(int j=0;j<m_numPeaks;j++){
		 out<<"set label \"p_{"<<j+1<<"}\" at ";
         for(int i=0;i<m_numDim;i++){
			 if(i+1<m_numDim)
            out<<mpp_peak[j][i]+0.1<<", ";
			 else out<<mpp_peak[j][i];
        }
        out<<endl;
    }*/

}

int DynamicContinuous::getNumofVisablePeaks(){
    return m_numVisablePeaks;

}
void DynamicContinuous::computeNumVisablePeaks(){
    m_numVisablePeaks=m_numPeaks;
    for(int i=0;i<m_numPeaks;i++){
		CodeVReal s(m_numDim,m_numObj);
		copy(mpp_peak[i],mpp_peak[i]+m_numDim,s.m_x.begin());
		evaluate_(s,false);
		double height=s.m_obj[0];
        switch(m_OptMode[0]){
            case MIN_OPT:
                if(height<mp_height[i]) m_numVisablePeaks--;
                break;
            case MAX_OPT:
                if(height>mp_height[i]) m_numVisablePeaks--;
                break;
        }
    }

}
bool DynamicContinuous::isVisable(const int rIdx){
	CodeVReal s(m_numDim,m_numObj);
	copy(mpp_peak[rIdx],mpp_peak[rIdx]+m_numDim,s.m_x.begin());
	evaluate_(s,false);
	double height=s.m_obj[0];
    switch(m_OptMode[0]){
        case MIN_OPT:
            if(height<mp_height[rIdx]) return false;
            break;
        case MAX_OPT:
            if(height>mp_height[rIdx]) return false;
            break;
    }
    return true;

}
void DynamicContinuous::addNoise(double *x_){
	for(int d=0; d<m_numDim;d++){
		double x=x_[d];
		x+=m_noiseSeverity_*Global::msp_global->mp_normalPro->Next();
		if(m_searchRange[d].m_upper<x) x=m_searchRange[d].m_upper;
		if(m_searchRange[d].m_lower>x)  x=m_searchRange[d].m_lower;
		x_[d]=x;
	}
}

int DynamicContinuous::getTrackNumber(int idex){
	return mp_isTracked[idex];
}

bool DynamicContinuous::isTracked(vector<double> &gen,vector<double> &obj){
	bool flag=false,movepeaks=false;
	for(int i=0;i<m_numPeaks;i++){
		double dis=0,dis1=fabs(obj[0]-mp_height[i]);
    	for(int j=0;j<m_numDim;j++) dis+=(gen[j]-mpp_peak[i][j])*(gen[j]-mpp_peak[i][j]);
		dis=sqrt(dis);
		if(dis<=m_disAccuracy&&dis1<=m_accuracy){
			// peak[i] assumed to be found
			int j=0;
			while(mp_heightOrder[j++]!=i&&j<m_numPeaks);
			if(!mp_found[i]){
				mp_isTracked[j-1]++;
				mp_found[i]=true;
				m_peaksFound++;
				flag=true;
			}
			
		}
		if(dis<m_timeLinkageSeverity){
			// move peak[i] to a near random position when it was tracked
			if(m_timeLinkageFlag){
				movePeak(i);
				updateTimeLinkage();
				movepeaks=true;
				m_flagTriggerTimeLinkage=true;
			}
		}
	}
	if(movepeaks){
		#ifdef DEMON_OFEC
		calculateSamplePoints();
		#endif
	}
	return flag;
}
bool DynamicContinuous::isTracked(double *gen,vector<double> &obj){
	bool flag=false,movepeaks=false;
	for(int i=0;i<m_numPeaks;i++){
		double dis=0,dis1=fabs(obj[0]-mp_height[i]);
    	for(int j=0;j<m_numDim;j++) dis+=(gen[j]-mpp_peak[i][j])*(gen[j]-mpp_peak[i][j]);
		dis=sqrt(dis);
		if(dis<=m_disAccuracy&&dis1<=m_accuracy){
			// peak[i] assumed to be found
			int j=0;
			while(mp_heightOrder[j++]!=i&&j<m_numPeaks);
			if(!mp_found[i]){
				mp_isTracked[j-1]++;
				mp_found[i]=true;
				m_peaksFound++;
				flag=true;
			}	
		}
		if(dis<m_timeLinkageSeverity){
			// move peak[i] to a near random position when it was tracked
			if(m_timeLinkageFlag){
				movePeak(i);
				updateTimeLinkage();
				movepeaks=true;
				m_flagTriggerTimeLinkage=true;
			}
		}
	}
	if(movepeaks){
		#ifdef DEMON_OFEC
		calculateSamplePoints();
		#endif
	}
	return flag;
}
int DynamicContinuous::getPeaksFound(){
	return m_peaksFound;
}
void DynamicContinuous::updateTimeLinkage(){
	if(!m_timeLinkageFlag) return;
	double range;
	for(int j=0;j<m_numDim;j++){
		range=fabs(double(m_searchRange[j].m_upper-m_searchRange[j].m_lower));
		mp_timeLinkage[j]=0.2*range*(Global::msp_global->mp_uniformPro->Next()-0.5);
	} 
}
void DynamicContinuous::movePeak(const int idx){
	if(idx<0||idx>=m_numPeaks) throw myException("index out of boundary @ DynamicContinuous::movePeak(const int idx)");
	for(int d=0; d<m_numDim;d++){
		double x=mpp_peak[idx][d];
		x+=mp_timeLinkage[d];
		if(m_searchRange[d].m_upper<x) x=m_searchRange[d].m_upper;
		if(m_searchRange[d].m_lower>x)  x=m_searchRange[d].m_lower;
		mpp_peak[idx][d]=x;
	}
}

double DynamicContinuous::getAssociateRadius(int idx){
	return mp_associateRadius[idx];
}
double DynamicContinuous::getPeaksTracedQaulity(){
	return m_peakQaulity;
}

bool DynamicContinuous::isGOptTracked(){
// the global optimum is assumed to be tracked if any one of the global optima is tracked	
	for(int i=0;i<m_numPeaks;i++){
		if(mp_globalOptimaIdx[i]&&mp_found[i]) return true;
	}
	return false;	
}
		//added 04/07/2014
bool DynamicContinuous::isGOpt(int idx){
// is peak i the global optimum	
	return mp_globalOptimaIdx[idx];
}

const double * DynamicContinuous::getNearestPeak(const vector<double>& p){
	int nearest=0;
	MyVector peak(GET_NUM_DIM,mpp_peak[0]);
	double dis=peak.getDis(p);
	for(int i=1;i<m_numPeaks;i++){
		copy(mpp_peak[i],mpp_peak[i]+GET_NUM_DIM,peak.begin());
		double d=peak.getDis(p);
		if(d<dis){
			dis=d;
			nearest=i;
		}
	}
	return mpp_peak[nearest];
}
