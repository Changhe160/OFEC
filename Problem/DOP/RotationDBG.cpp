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

#include "RotationDBG.h"
#include "../../Global/global.h"

RotationDBG::RotationDBG(ParamMap &v):Problem((v[param_proId]),(v[param_numDim]),string(),1),\
	RealDBG((v[param_proId]),(v[param_numDim]),\
	(v[param_numPeak])){
	 m_name="DYN_CONT_RotationDBG";
	 addProTag(MMP);
	initialize(static_cast<ChangeType>((int)(v[param_changeType])),(v[param_changeRatio]),(v[param_flagNumDimChange]),\
		(v[param_flagNumPeakChange]),\
		(v[param_peakNumChangeMode]),(v[param_flagNoise]),(v[param_flagTimeLinkage]));
}
RotationDBG::RotationDBG(const int rId, const int rDimNumber, const int rNumPeaks,\
                         const ChangeType rT, const double rChangingRatio, const bool rFlagDimChange,const bool rFlagNumPeakChange,\
						 const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage):\
						 Problem(rId,rDimNumber,string(),1),RealDBG(rId,rDimNumber,rNumPeaks){
    //ctor
    m_name="DYN_CONT_RotationDBG";
	addProTag(MMP);
    initialize(rT,rChangingRatio,rFlagDimChange,rFlagNumPeakChange,peakNumChangeMode,flagNoise, flagTimelinkage);
}

RotationDBG::~RotationDBG()
{
    //dtor
}

void  RotationDBG::setWidth(const double w){

    for(int i=0;i<m_numPeaks;i++)
        if(m_changeType.type==CT_Chaotic)
            mp_width[i]=m_minWidth+(m_maxWidth-m_minWidth)*Global::msp_global->mp_uniformPro->Next();
        else
            mp_width[i]=w;
};

RotationDBG& RotationDBG::operator =(const RotationDBG & r_dbg){
	if(this==& r_dbg) return *this;
	if(m_numDim!=r_dbg.m_numDim||m_numPeaks!=r_dbg.m_numPeaks) throw myException("RotationDBG  assignment@RotationDBG::operator =");
	RealDBG::operator =(r_dbg);
	return *this;
}
void RotationDBG::randomChange(){
	for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	heightStandardChange();
	widthStandardChange();
	positionStandardChange(0);

	restoreInfor();
    calculateGlobalOptima();
     updateNumberofChanges();
	m_changeType.counter++;
}
void RotationDBG::recurrentChange(){
    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);


	double initial_angle;

	double height_range=m_maxHeight-m_minHeight;
	double width_range=m_maxWidth-m_minWidth;
	for(int i=0;i<m_numPeaks;i++){
	    if(mp_whetherChange[i]==false) continue;
		initial_angle=(double)getPeriod()*i/m_numPeaks;
		mp_height[i]=m_minHeight+height_range*(sin(2*OFEC_PI*(m_changeType.counter+initial_angle)/getPeriod())+1)/2.;
		mp_width[i]=m_minWidth+width_range*(sin(2*OFEC_PI*(m_changeType.counter+initial_angle)/getPeriod())+1)/2.;
	}
	initial_angle=OFEC_PI*(sin(2*OFEC_PI*(m_changeType.counter)/getPeriod())+1)/12.;
	positionStandardChange(initial_angle);

    restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void RotationDBG::chaoticChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	for(int i=0;i<m_numPeaks;i++){
	    if(mp_whetherChange[i]==false) continue;
		mp_height[i]=gChaoticValue(mp_height[i],m_minHeight,m_maxHeight);
		mp_width[i]=gChaoticValue(mp_width[i],m_minWidth,m_maxWidth);
	}
	positionStandardChange(0);

	restoreInfor();
    calculateGlobalOptima();
     updateNumberofChanges();
	m_changeType.counter++;
}
void RotationDBG::smallStepChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	heightStandardChange();
	widthStandardChange();
	positionStandardChange(0);

	restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void RotationDBG::largeStepChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	heightStandardChange();
	widthStandardChange();
	positionStandardChange(0);

	restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void RotationDBG::recurrentNoisyChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	double initial_angle;
	double height_range=m_maxHeight-m_minHeight;
	double width_range=m_maxWidth-m_minWidth;
	double noisy;
	for(int i=0;i<m_numPeaks;i++){
        if(mp_whetherChange[i]==false) continue;
		initial_angle=(double)getPeriod()*i/m_numPeaks;
		mp_height[i]=sinValueNoisy(m_changeType.counter,m_minHeight,m_maxHeight,height_range,initial_angle,m_noisySeverity);
		mp_width[i]=sinValueNoisy(m_changeType.counter,m_minWidth,m_maxWidth,width_range,initial_angle,m_noisySeverity);
	}
	initial_angle=OFEC_PI*(sin(2*OFEC_PI*(m_changeType.counter)/getPeriod())+1)/12.;

	noisy=m_noisySeverity*Global::msp_global->mp_normalPro->Next();
	positionStandardChange(initial_angle+noisy);

    restoreInfor();
	calculateGlobalOptima();
    updateNumberofChanges();
	m_changeType.counter++;
}
ReturnFlag  RotationDBG::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag2){
	CodeVReal &s=dynamic_cast<CodeVReal &>(ss);
	double *gene=new double[m_numDim];
	copy(s.m_x.begin(),s.m_x.end(),gene);

	if(this->m_noiseFlag)	addNoise(gene);
	vector<double> fit(m_numPeaks,0);
	for(int i=0;i<m_numPeaks;i++){
		for(int j=0;j<m_numDim;j++)
			fit[i]+=(gene[j]-mpp_peak[i][j])*(gene[j]-mpp_peak[i][j]);
		if(fit[i]!=0) fit[i]=sqrt(fit[i]/m_numDim);
		fit[i]=mp_height[i]/(1+mp_width[i]*fit[i]);
	}
	s.m_obj[0]=*max_element(fit.begin(),fit.end());
	if(rFlag&&m_evals%m_changeFre==0) Solution<CodeVReal>::initilizeWB(s);
	if(rFlag){
		isTracked(gene,s.m_obj);
		m_evals++;
	}
	bool flag;
	#ifdef OFEC_CONSOLE
		if(Global::msp_global->mp_algorithm!=nullptr)	flag=!Global::msp_global->mp_algorithm->ifTerminating();
		else flag=true;
	#endif
	#ifdef OFEC_DEMON
		flag=true;
	#endif
    if(rFlag&&m_evals%m_changeFre==0&&flag){
		DynamicProblem::change();
	} 
	delete []gene;
	gene=0;
	ReturnFlag rf=Return_Normal;
	if(rFlag){
		if(Global::msp_global->mp_algorithm!=nullptr){
				if(Global::msp_global->mp_algorithm->ifTerminating()){ rf=Return_Terminate; }
				else if(Global::msp_global->mp_problem->isProTag(DOP)){	
					if(CAST_PROBLEM_DYN->getFlagTimeLinkage()&&CAST_PROBLEM_DYN->getTriggerTimelinkage()){
						rf=Return_Change_Timelinkage; 
					}
					if((Global::msp_global->mp_problem->getEvaluations()+1)%(CAST_PROBLEM_DYN->getChangeFre())==0){
						rf=Return_ChangeNextEval; 
					}
					if(Global::msp_global->mp_problem->getEvaluations()%(CAST_PROBLEM_DYN->getChangeFre())==0){
						if(CAST_PROBLEM_DYN->getFlagDimensionChange()){
							rf=Return_Change_Dim; 
						}
						rf=Return_Change;
					}
				}
		}
	}
	return rf;
}


void RotationDBG::widthStandardChange(){
	double step;
	for(int i=0;i<m_numPeaks;i++){
        if(mp_whetherChange[i]==false) continue;
		step=mp_widthSeverity[i]*standardChange(m_changeType.type,m_minWidth,m_maxWidth);
		mp_width[i]=mp_width[i]+step;

		if(mp_width[i]>m_maxWidth||mp_width[i]<m_minWidth) mp_width[i]=mp_width[i]-step;

	}
}

void RotationDBG::changeDimension(){
    /// no need to preserve  previous information, e.g., positions, heights, width....


	RotationDBG* r_dbg=new RotationDBG(Global::msm_pro["DYN_CONT_RotationDBG"],m_dimNumberTemp,m_numPeaks,m_changeType.type,
                                    m_changePeakRatio,m_flagDimensionChange,m_flagNumPeaksChange,
									m_numPeaksChangeMode,m_noiseFlag,m_timeLinkageFlag);
	if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy)
		r_dbg->setPeriod(m_period);

	r_dbg->parameterSetting(this);
	r_dbg->calculateGlobalOptima();

    RealDBG::freeMemory();
	DynamicContinuous::freeMemory();
    Problem::freeMemory();

	Problem::allocateMemory(m_dimNumberTemp);
	ContinuousProblem::allocateMemory(m_dimNumberTemp);
    DynamicContinuous::allocateMemory(m_dimNumberTemp,m_numPeaks);
    RealDBG::allocateMemory(m_dimNumberTemp,m_numPeaks);
	
    m_numDim=m_dimNumberTemp;
    setPeriod(m_period);
	*this=*r_dbg;
	delete r_dbg;
	r_dbg=0;
}

void RotationDBG::parameterSetting(Problem * rP){
	RealDBG::parameterSetting(rP);
}

void RotationDBG::initialize(const ChangeType rT,double const rChangingRatio, const bool rFlagDimChange,const bool rFlagNumPeakChange,
	const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage){
    RealDBG::initialize(rT,rFlagDimChange,rFlagNumPeakChange,peakNumChangeMode,flagNoise,flagTimelinkage);
    
	m_globalOpt.setFlagLocTrue();
    setOptType(MAX_OPT);
	calculateGlobalOptima();

	setNumberofChanges(rChangingRatio);
	setAccuracy(0.1);
	setDisAccuracy(0.1);
	
}


void  RotationDBG::changeNumPeaks(){

    RotationDBG* r_dbg=new RotationDBG(Global::msm_pro["DYN_CONT_RotationDBG"],m_numDim,m_numPeaksTemp,m_changeType.type,\
									m_changePeakRatio,m_flagDimensionChange,m_flagNumPeaksChange,\
									m_numPeaksChangeMode,m_noiseFlag,m_timeLinkageFlag);

	if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy)
		r_dbg->setPeriod(m_period);
	r_dbg->parameterSetting(this);
	r_dbg->calculateGlobalOptima();


    RealDBG::freeMemory();
	DynamicContinuous::freeMemory();

    DynamicContinuous::allocateMemory(m_numDim,m_numPeaksTemp);
    RealDBG::allocateMemory(m_numDim,m_numPeaksTemp);

    m_numPeaks=m_numPeaksTemp;
    setPeriod(m_period);
	*this=*r_dbg;
	delete r_dbg;

}

