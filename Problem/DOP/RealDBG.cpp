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

#include "RealDBG.h"
#include "../../Global/global.h"

RealDBG::RealDBG(const int rId, const int rDimNumber,  const int rNumPeaks, const int numObj):DynamicContinuous(rId,rDimNumber,rNumPeaks,numObj)
{
    //ctor
    allocateMemory(m_numDim,m_numPeaks);
    mppp_rotationPlanes=0;
	setSearchRange(-5,5);
   
    m_maxHeight=100;
    m_minHeight=10;
    m_maxWidth=10;
    m_minWidth=1;
    setNoisySeverity(0.8);
    setHeightSeverity(5.0);
    setWidthSeverity(0.5);
    setWidth(5); /// value in 1-10


}

void RealDBG::allocateMemory(const int rDimNum, const int rPeaks){
     mp_rotationMatrix= new Matrix [rPeaks];
	 for(auto i=0;i<rPeaks;i++){
		 mp_rotationMatrix[i].resize(rDimNum,rDimNum);
	 }
}

RealDBG::~RealDBG(){
    //dtor
    freeMemory();
}

RealDBG & RealDBG::operator=(const RealDBG & rP){

    if(this==&rP) return *this;

    DynamicContinuous::operator=(rP);
    copy(rP.mp_rotationMatrix,rP.mp_rotationMatrix+m_numPeaks,mp_rotationMatrix);

    m_prediction=rP.m_prediction;

    if(getPeriod()!=rP.getPeriod()){
        throw myException("The period must be the same@RealDBG::operator=");
    }

    for(int i=0;i<m_period;i++){
		for(int j=0;j<m_numPeaks;j++)
            copy(rP.mppp_rotationPlanes[i][j],rP.mppp_rotationPlanes[i][j]+m_numDim,mppp_rotationPlanes[i][j]);
	}

    return *this;

}
void  RealDBG::freeMemory(){
    delete[] mp_rotationMatrix;

    if(mppp_rotationPlanes){
		for(int j=0;j<getPeriod();j++){
			for(int i=0;i<m_numPeaks;i++)
				delete [] mppp_rotationPlanes[j][i];
			delete []mppp_rotationPlanes[j];
		}
		delete []mppp_rotationPlanes;
	}

	mp_rotationMatrix=0;
	mppp_rotationPlanes=0;
}
bool RealDBG::setPeriod(const int p){
	if(p<1) return false;
	DynamicProblem::setPeriod(p);

	mppp_rotationPlanes=new int**[m_period];
	for(int i=0;i<m_period;i++){
		mppp_rotationPlanes[i]=new int*[m_numPeaks];
		for(int j=0;j<m_numPeaks;j++)
			mppp_rotationPlanes[i][j]=new int[m_numDim];
	}
	return true;
}

void RealDBG::correctSolution(double *x){
	for(int j=0;j<m_numDim;j++){
		if(m_searchRange[j].m_upper<x[j])
			x[j]=m_searchRange[j].m_upper;
		else if(m_searchRange[j].m_lower>x[j])
			x[j]=m_searchRange[j].m_lower;
	}
}
void RealDBG::heightStandardChange(){

	double step;
	for(int i=0;i<m_numPeaks;i++){
		if(mp_whetherChange[i]==false) continue;
		step=mp_heightSeverity[i]*standardChange(getChangeType(),m_minHeight,m_maxHeight);
		mp_height[i]=mp_height[i]+step;
		if(mp_height[i]>m_maxHeight||mp_height[i]<m_minHeight) mp_height[i]=mp_height[i]-step;

	}
}
void RealDBG::positionStandardChange(double angle){

	
	if(getChangeType()==CT_Chaotic){
		for(int i=0;i<m_numPeaks;i++){
		    if(mp_whetherChange[i]==false) continue;
			for(int j=0;j<m_numDim;j++)
				mpp_peak[i][j]=gChaoticValue(mpp_peak[i][j],m_searchRange[j].m_lower,m_searchRange[j].m_upper,m_chaoticConstant);
		}
		return;
	}
	
	// for each basic function of dimension n(even number) , R=R(l1,l2)*R(l3,l4)*....*R(li-1,li), 0<=li<=n

	int * d=new int[m_numDim];
	Matrix I(m_numDim,m_numDim);
	double *gene=new double[m_numDim];
	for(int i=0;i<m_numPeaks;i++){
		if((getChangeType()==CT_Recurrent||getChangeType()==CT_RecurrentNoisy)&&m_changeType.counter>=getPeriod()){
			copy(mppp_rotationPlanes[m_changeType.counter%getPeriod()][i],mppp_rotationPlanes[m_changeType.counter%getPeriod()][i]+m_numDim,d);
		}
		else{
			Global::msp_global->initializeRandomArray(d,m_numDim,Program_Problem);
			if(getChangeType()==CT_Recurrent||getChangeType()==CT_RecurrentNoisy)
			copy(d,d+m_numDim,mppp_rotationPlanes[m_changeType.counter][i]);
		}

		if((getChangeType()==CT_Recurrent||getChangeType()==CT_RecurrentNoisy)&&m_changeType.counter%getPeriod()==0)
			copy(mpp_initialPeak[i],mpp_initialPeak[i]+m_numDim,mpp_peak[i]);

		if(mp_whetherChange[i]==false) continue;

		I.identity();
		for(int j=0;j+1<m_numDim;j+=2){
			if(getChangeType()==CT_SmallStep||getChangeType()==CT_LargeStep||getChangeType()==CT_Random)
				angle=standardChange(getChangeType(), -OFEC_PI,OFEC_PI);
			I.setRotationAngle(d[j],d[j+1],angle);
			if(j==0) mp_rotationMatrix[i]=I;
			else
				mp_rotationMatrix[i]*=I;
		}
		Matrix m(m_numDim,1);
		m.setDataRow(mpp_peak[i],m_numDim);
		m*=mp_rotationMatrix[i];
		copy(m[0].begin(),m[0].end(),gene);
		correctSolution(gene);
		copy(gene,gene+m_numDim,mpp_peak[i]);
	}
	delete [] d;
	delete []gene;
	d=0;
	gene=0;
}


void RealDBG::parameterSetting(Problem * rP){

    DynamicContinuous::parameterSetting(rP);

    RealDBG *r_dbg=dynamic_cast<RealDBG *>(rP);
	int dim=m_dimNumberTemp<rP->getNumDim()?m_dimNumberTemp:rP->getNumDim();
   int peaks=m_numPeaks<r_dbg->getNumberofPeak()?m_numPeaks:r_dbg->getNumberofPeak();

    m_prediction=r_dbg->m_prediction;

    if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy){
		for(int i=0;i<r_dbg->m_period;i++){
			if(m_changeType.counter<=i) break;
			for(int j=0;j<peaks;j++){
				if(dim==m_dimNumberTemp){// the number of dimensions decreases
					for(int m=0,k=0;k<dim;k++,m++)
						if(r_dbg->mppp_rotationPlanes[i][j][m]==dim) {k--;continue;}
						else
							mppp_rotationPlanes[i][j][k]=r_dbg->mppp_rotationPlanes[i][j][m];

				}else
					copy(r_dbg->mppp_rotationPlanes[i][j],r_dbg->mppp_rotationPlanes[i][j]+dim,mppp_rotationPlanes[i][j]);

			}

		}
	}
}

double  RealDBG::standardChange(const ChangeType T, const double min, const double max){
	double step,sign;
	switch(T){
		case CT_SmallStep:
			step=-1+2*Global::msp_global->mp_uniformPro->Next();
			step=m_alpha*step*(max-min);
			break;
		case CT_Random:
			step=Global::msp_global->mp_normalPro->Next();
			break;
		case CT_LargeStep:
			step=-1+2*Global::msp_global->mp_uniformPro->Next();
			if(step>0)sign=1;
			else if(step<0) sign=-1;
			else sign=0;
			step=(m_alpha*sign+(m_maxAlpha-m_alpha)*step)*(max-min);
			break;
		case CT_Recurrent:
		case CT_Chaotic:
		case CT_RecurrentNoisy:
			break;
		}
	return step;
}
void RealDBG::reinitialize(){
	if(m_numPeaks!= *DynamicProblem::ms_initNumPeaks) { 
		m_numPeaksTemp= *DynamicProblem::ms_initNumPeaks;  
		changeNumPeaks();
	}
	if(m_numDim!=*DynamicProblem::ms_initNumDim){ 
		m_dimNumberTemp=*DynamicProblem::ms_initNumDim; 
		changeDimension();
	}

}
void RealDBG::reset(){

    m_changeType.counter=0;
    m_changeCounter=0;
    double *t=new double[m_numPeaks];

	for(int i=0;i<m_numPeaks;i++){
		if(m_changeType.type==CT_Chaotic)
			t[i]=m_minHeight+(m_maxHeight-m_minHeight)*Global::msp_global->mp_uniformPro->Next();
		else
			t[i]=50;
	}
	setHeight(t);
    delete [] t;


	double **position;
    position=new double*[m_numPeaks];
    for(int i=0;i<m_numPeaks;i++)
        position[i]=new double[m_numDim];
    for(int i=0;i<m_numPeaks;i++){
        for(int j=0;j<m_numDim;j++){
            position[i][j]=m_searchRange[j].m_lower+(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*Global::msp_global->mp_uniformPro->Next();
        }
    }
    setPosition(position);


	for(int i=0;i<m_numPeaks;i++)		delete []position[i];
	delete [] position;
	position=0;

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	calculateGlobalOptima();
}
void RealDBG::initialize(const ChangeType rT, const bool rFlagDimChange, const bool rFlagNumPeakChange
	,const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage){
   
    setDimensionChange(rFlagDimChange);
	setNumPeakChangeMode(peakNumChangeMode);
	setNumPeaksChange(rFlagNumPeakChange);
    setNoiseFlag(flagNoise);
	setTimeLinkageFlag(flagTimelinkage);


    if(m_flagDimensionChange){
        m_changeType.type=CT_Random;
        m_dirDimensionChange=true;
    }else if(m_flagNumPeaksChange){
         m_changeType.type=CT_Random;
         m_dirNumPeaksChange=true;
	}else if(m_noiseFlag||m_timeLinkageFlag){
		 m_changeType.type=CT_Random;
         
	}else{
         m_changeType.type=rT;
    }
    m_changeType.counter=0;

    if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy)      setPeriod(12);
    else      setPeriod(0);


    double *t=new double[m_numPeaks];

    setChoaticConstant(3.67);

	for(int i=0;i<m_numPeaks;i++){
		if(m_changeType.type==CT_Chaotic)
			t[i]=m_minHeight+(m_maxHeight-m_minHeight)*Global::msp_global->mp_uniformPro->Next();
		else
			t[i]=50;
	}
	setHeight(t);
    delete [] t;
	t=0;

	double **position;
    position=new double*[m_numPeaks];
    for(int i=0;i<m_numPeaks;i++)
        position[i]=new double[m_numDim];
    for(int i=0;i<m_numPeaks;i++){
        for(int j=0;j<m_numDim;j++){
            position[i][j]=m_searchRange[j].m_lower+(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*Global::msp_global->mp_uniformPro->Next();
        }
    }
    setPosition(position);


	for(int i=0;i<m_numPeaks;i++)		delete []position[i];
	delete [] position;
	position=0;
	for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

}

 void RealDBG::restoreInfor(){

     for(int i=0;i<m_numPeaks;i++){
        if(!mp_whetherChange[i]){
            copy(mpp_prePeak[i],mpp_prePeak[i]+m_numDim,mpp_peak[i]);
            mp_height[i]=mp_preHeight[i];
            mp_width[i]=mp_preWidth[i];
        }
     }

  }
