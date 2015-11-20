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

#include "DynamicProblem.h"
#include "../../Global/global.h"
#include "../../Measure/mSingleObj.h"

boost::thread_specific_ptr<int> DynamicProblem::ms_initNumPeaks,DynamicProblem::ms_initNumDim,DynamicProblem::ms_numInstance;

#ifdef OFEC_DEMON
#include "../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

DynamicProblem::DynamicProblem(const int rId, const int rDimNumber, const int rNumPeaks,const int runId,const unsigned numObj):Problem(rId,rDimNumber, string(),numObj),m_changeCounter(0)
,m_dimNumberTemp(rDimNumber),m_numPeaks(rNumPeaks),m_numPeaksTemp(rNumPeaks),m_noiseFlag(false),m_timeLinkageFlag(false),m_flagTriggerTimeLinkage(false){
    //ctor

    m_changeFre=5000;
    m_changeType.type=CT_Random;
    m_changeType.counter=0;
    m_period=0;
    m_flagDimensionChange=false;
    m_dirDimensionChange=true;
    m_synchronize=true;
    m_noisySeverity=0.8;

    m_alpha=0.04;
    m_maxAlpha=0.1;
    m_chaoticConstant=3.67; //in [3.57,4]

    m_flagNumPeaksChange=false;
    m_dirNumPeaksChange=true;
	m_numPeaksChangeMode=1;
	m_noiseSeverity_=0.01;
	m_timeLinkageSeverity=0.1;
	m_proPar<<"Change frequency:"<<m_changeFre<<"; "<<"TotalEvals:"<<Global::g_arg[param_maxEvals]<<"; "<<"Peaks:"<<m_numPeaks<<"; "<<"NumPeaksChange:"<<m_flagNumPeaksChange<<"-"<<m_numPeaksChangeMode<<"; "<<
		"NoisyEnvioronments:"<<m_noiseFlag<<"; NoiseSeverity:"<<m_noiseSeverity_<<"; TimeLinkageEnvironments:"<<m_timeLinkageFlag<<"; TimeLinkageSeverity:"<<m_timeLinkageSeverity<<"; DimensionalChange:"<<m_flagDimensionChange<<"; ";

	if(!ms_numInstance.get()) ms_numInstance.reset(new int(0));
	if(!ms_initNumPeaks.get()) ms_initNumPeaks.reset(new int);
	if(!ms_initNumDim.get())ms_initNumDim.reset(new int);
	(*ms_numInstance)++;
	if(*ms_numInstance==1){
		*ms_initNumPeaks=m_numPeaks;
		*ms_initNumDim=m_numDim;
	}
	addProTag(DOP);
}

DynamicProblem::~DynamicProblem()
{
    //dtor
}
int DynamicProblem::getNumberofPeak()const{
	return m_numPeaks;
}

void DynamicProblem::setNumPeaksChange(const bool rPC){
	m_flagNumPeaksChange=rPC;
    size_t start, end;
    start=m_proPar.str().find("NumPeaksChange:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
	
    stringstream ss;
    ss<<"NumPeaksChange:"<<m_flagNumPeaksChange<<"-"<<m_numPeaksChangeMode<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);
}

void DynamicProblem::setChangeFre(const int rChangeFre){
    if(rChangeFre>0) m_changeFre=rChangeFre;
    else{
        throw myException("Change frequncy must be greater than 0 @DynamicProblem::setChangeFre");
    }

    size_t start, end;
    start=m_proPar.str().find("Change frequency:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Change frequency:"<<m_changeFre<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);


}
bool DynamicProblem::setPeriod(const int rPeriod){
    if(rPeriod>=0) m_period=rPeriod;
    else{
        throw myException("period must be positive@ DynamicProblem::setPeriod");
    }
    return true;
}
void DynamicProblem::setChangeType(const SChangeType &rChangeType){
    m_changeType=rChangeType;
}
void DynamicProblem::setChangeType(const ChangeType rT){
    m_changeType.type=rT;
}

void DynamicProblem::setDimensionChange(const bool rFlag){

    m_flagDimensionChange=rFlag;

	size_t start, end;
    start=m_proPar.str().find("DimensionalChange:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"DimensionalChange:"<<m_flagDimensionChange<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);

}
void DynamicProblem::setNoiseSeverity_(double value){ 
	m_noiseSeverity_=value;
	size_t start, end;
    start=m_proPar.str().find("NoiseSeverity:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"NoiseSeverity:"<<m_noiseSeverity_<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);
}
void DynamicProblem::setTimeLinkageSeverity(double value){
	m_timeLinkageSeverity=value;
	
	size_t start, end;
    start=m_proPar.str().find("TimeLinkageSeverity:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"TimeLinkageSeverity:"<<m_timeLinkageSeverity<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);
}


void DynamicProblem::setChangeDirction(const bool rFlag){
    m_dirDimensionChange=rFlag;
}
void DynamicProblem::setSynchronize(const bool rFlag){
    m_synchronize=rFlag;
}

void DynamicProblem::setNoisySeverity(const double rSeverity){

    m_noisySeverity=rSeverity;
}


void DynamicProblem::change(){

	m_changeCounter++;
	switch(getChangeType()){
	case CT_Random:
		randomChange();
		break;
	case CT_Recurrent:
		recurrentChange();
		break;
	case CT_RecurrentNoisy:
		recurrentNoisyChange();
		break;
	case CT_SmallStep:
		smallStepChange();
		break;
	case CT_LargeStep:
		largeStepChange();
		break;
	case CT_Chaotic:
		chaoticChange();
		break;
	default :
		break;
	}

    if(m_flagDimensionChange){
		
        if(m_numDim==msc_MinDimensionNumber)
        m_dirDimensionChange=true;
        if(m_numDim==msc_MaxDimensionNumber)
        m_dirDimensionChange=false;

        if(m_dirDimensionChange==true) {
			m_dimNumberTemp+=1;
        }else {
			m_dimNumberTemp-=1;
        }
        changeDimension();
    }

    if(m_flagNumPeaksChange){
		if(m_numPeaksChangeMode==1||m_numPeaksChangeMode==2){
			if((unsigned int)m_numPeaks>=msc_MaxNumPeaks-1) m_dirNumPeaksChange=false;
			if((unsigned int)m_numPeaks<=msc_MinNumPeaks+1) m_dirNumPeaksChange=true;
			int step=0;
			 
			if(IS_PROBLEM_NAME(m_id,"DYN_CONT_CompositionDBG")) step=2;
			else if(IS_PROBLEM_NAME(m_id,"DYN_CONT_RotationDBG")) step=2;
			else if(IS_PROBLEM_NAME(m_id,"DYN_CONT_MovingPeak")) step=2;
			else step=2;

			if(m_numPeaksChangeMode==2){
				step=Global::msp_global->getRandInt(step/2,5*step/2,Program_Problem);
			}
			
			if(m_dirNumPeaksChange==true){
				if(m_numPeaks+step<=msc_MaxNumPeaks)		m_numPeaksTemp=m_numPeaks+step;
				else m_numPeaksTemp=msc_MaxNumPeaks;
			}else{ 
				if(m_numPeaks-step>=msc_MinNumPeaks)		m_numPeaksTemp=m_numPeaks-step;
				else m_numPeaksTemp=msc_MinNumPeaks;
			}
		}else{
			//random change
			m_numPeaksTemp=Global::msp_global->getRandInt(msc_MinNumPeaks,msc_MaxNumPeaks,Program_Problem);
		}
        changeNumPeaks();
    }

	
    if(mSingleObj::getSingleObj()!=nullptr){
		vector<double> gOpt;
		if(Global::msp_global->mp_problem->getObjGlobalOpt(gOpt)){
		   mSingleObj::getSingleObj()->addGOpt(Global::msp_global->m_runId,gOpt[0]);
		}else{
			cout<<"err"<<endl;
		}
	}
	
	#ifdef OFEC_DEMON
		msp_buffer->updateFitnessLandsacpe_();
	#endif
	
}



 DynamicProblem & DynamicProblem::operator=(const DynamicProblem & rDP){
    if(this==&rDP) return *this;

    if(m_numDim!=rDP.m_numDim){
        throw myException("The number of dimensions must be same!@DynamicProblem::operator=");
    }
    if(m_changeType.type!=rDP.m_changeType.type){
        throw myException("The change type must be same!@DynamicProblem::operator=");
    }
    if(m_numPeaks!=rDP.m_numPeaks){
        throw myException("The number of peaks must be same!@DynamicProblem::operator=");
    }
    Problem::operator=(rDP);


    m_changeType.counter=rDP.m_changeType.counter;
    m_changeFre=rDP.m_changeFre;
    m_period=rDP.m_period;
    m_flagDimensionChange=rDP.m_flagDimensionChange;
    m_dirDimensionChange=rDP.m_dirDimensionChange;
    m_synchronize=rDP.m_synchronize;
    m_noisySeverity=rDP.m_noisySeverity;
   
    m_alpha =rDP.m_alpha;
    m_maxAlpha= rDP.m_maxAlpha;
    m_chaoticConstant =rDP.m_chaoticConstant;

    m_flagNumPeaksChange=rDP.m_flagNumPeaksChange;
    m_dirNumPeaksChange=rDP.m_dirNumPeaksChange;
	m_numPeaksChangeMode=rDP.m_numPeaksChangeMode;

	m_noiseFlag=rDP.m_noiseFlag;
	m_timeLinkageFlag=rDP.m_timeLinkageFlag;

	m_noiseSeverity_=rDP.m_noiseSeverity_;
	m_timeLinkageSeverity=rDP.m_timeLinkageSeverity;
		m_flagTriggerTimeLinkage=rDP.m_flagTriggerTimeLinkage;
    return *this;
 }

 void  DynamicProblem::parameterSetting(Problem * rdp){

	 Problem::parameterSetting(rdp);

	 DynamicProblem *rDP=dynamic_cast<DynamicProblem*>(rdp);
    m_changeType=rDP->m_changeType;
    m_changeFre=rDP->m_changeFre;
    m_period=rDP->m_period;
    m_flagDimensionChange=rDP->m_flagDimensionChange;
    m_dirDimensionChange=rDP->m_dirDimensionChange;
    m_synchronize=rDP->m_synchronize;
    m_noisySeverity=rDP->m_noisySeverity;
    m_alpha =rDP->m_alpha;
    m_maxAlpha= rDP->m_maxAlpha;
    m_chaoticConstant =rDP->m_chaoticConstant;

    m_flagNumPeaksChange=rDP->m_flagNumPeaksChange;
    m_dirNumPeaksChange=rDP->m_dirNumPeaksChange;
	m_numPeaksChangeMode=rDP->m_numPeaksChangeMode;
	m_noiseFlag=rDP->m_noiseFlag;
	m_timeLinkageFlag=rDP->m_timeLinkageFlag;
	m_noiseSeverity_=rDP->m_noiseSeverity_;
	m_timeLinkageSeverity=rDP->m_timeLinkageSeverity;
		m_flagTriggerTimeLinkage=rDP->m_flagTriggerTimeLinkage;
 }

 double DynamicProblem::sinValueNoisy(const int x,const double min, const double max, const double amplitude, const double angle,const double noisy_severity){
														// return a value in recurrent with noisy dynamism environment
	double y;
	double noisy,t;
	y=min+amplitude*(sin(2*OFEC_PI*(x+angle)/m_period)+1)/2.;
	noisy=noisy_severity*Global::msp_global->mp_normalPro->Next();
	t=y+noisy;
	if(t>min&&t<max) y=t;
	else y= t-noisy;
	return y;
}

double DynamicProblem::chaoticStep(const double x, const double min, const double max, const double scale){
	if(min>max) return -1;
	double chaotic_value;
	chaotic_value=(x-min)/(max-min);
	chaotic_value=m_chaoticConstant*chaotic_value*(1-chaotic_value);
	//return fabs((min+chaotic_value*(max-min)-x)* Global::scale);
	return chaotic_value*scale;
}

bool DynamicProblem::predictChange(const int evalsMore){
    int fre=getChangeFre();
    int evals=getEvaluations()%fre;
    if(evals+evalsMore>=fre) return true;
    else return false;
}
void DynamicProblem::setNumPeakChangeMode(const int mode){
	m_numPeaksChangeMode=mode;
}
int DynamicProblem::getNumPeakChangeMode(){
	return m_numPeaksChangeMode;
}
void DynamicProblem::setNoiseFlag(const bool flag){
	m_noiseFlag=flag;

	size_t start, end;
    start=m_proPar.str().find("NoisyEnvioronments:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
	
    stringstream ss;
    ss<<"NoisyEnvioronments:"<<m_noiseFlag<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);

} 
void DynamicProblem::setTimeLinkageFlag(const bool flag){
	m_timeLinkageFlag=flag;
	size_t start, end;
    start=m_proPar.str().find("TimeLinkageEnvironments:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
	
    stringstream ss;
    ss<<"TimeLinkageEnvironments:"<<m_timeLinkageFlag<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);
}
int DynamicProblem::getInitialNumPeaks(){
	return *ms_initNumPeaks;
}

