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
#include "../../Global/global.h"
#include "MovingPeak.h"

static double basis_peak[5][7]={
  {8.0,  64.0,  67.0,  55.0,   4.0, 0.1, 50.0},
  {50.0,  13.0,  76.0,  15.0,   7.0, 0.1, 50.0},
  {9.0,  19.0,  27.0,  67.0, 24.0, 0.1, 50.0},
  {66.0,  87.0,  65.0,  19.0,  43.0, 0.1, 50.0},
  {76.0,  32.0,  43.0,  54.0,  65.0, 0.1, 50.0}
};

static double twin_peak [7] = /* difference to first peak */
  {
    1.0,  1.0,  1.0,  1.0,   1.0, 0.0, 0.0,
  };
void  MovingPeak::freeMemory(){
	int i;

	for (i=0; i< m_numPeaks; i++){
        delete [] mpp_prevMovement[i];
	}
	delete [] mpp_prevMovement;
	delete []  mp_shift ;
	delete []  mp_coveredPeaks ;

	mpp_prevMovement=0;
	mp_shift=0 ;
	mp_coveredPeaks=0 ;


}
void MovingPeak::allocateMemory(const int rDimNum, const int rPeaks){
    mp_shift = new double[rDimNum];
	mp_coveredPeaks = new int[rPeaks];
	mpp_prevMovement = new double*[rPeaks];

	for (int i=0; i< rPeaks; i++){
        mpp_prevMovement[i] =new double[rDimNum];
	}

}


bool MovingPeak::readData(){

	/***************************************
		//		m_F		Evaluation Function
		//		1		constant_basis_func()
		//		2		five_peak_basis_func()
		//		3		peak_function1()
		//		4		peak_function_cone()
		//		5		peak_function_hilly()
		//		6		peak_function_twin()
     **************************************
		in>>temp>>m_vlength; // distance by which the peaks are moved, severity
		 lambda determines whether there is a direction of the movement, or whether
		they are totally random. For lambda = 1.0 each move has the same direction,
		while for lambda = 0.0, each move has a random direction
		//in>>temp>>lambda;
		//in>>temp>>m_useBasisFunction;  if set to 1, a static landscape (basis_function) is included in the fitness evaluation
	}*/

	m_F=4;
	m_lambda=0;
	m_useBasisFunction=0;
	m_vlength=1.0;
	m_calculateRightPeak = 1; /* saves computation time if not needed and set to 0 */
	m_minHeight = 30.0;
	m_maxHeight = 70.0;
	m_standardHeight = 50.0;
	m_minWidth = 1.0 ;
	m_maxWidth = 12.0;
	m_standardWidth = 0.0;

	setSeverity();
	return true;

}

void MovingPeak::setSeverity(){
	for (int i=0; i< m_numPeaks; i++){
		mp_heightSeverity[i]= Global::msp_global->getRandFloat(1,10,Program_Problem);//1.+9.*i/m_numPeaks;//severity of height changes, larger numbers  mean larger severity. in the contex of ROOT, peaks have different values
		mp_widthSeverity[i]=Global::msp_global->getRandFloat(0.1,1,Program_Problem);// 0.1+0.9*i/m_numPeaks;severity of width changes, larger numbers mean larger severity
	}

}

MovingPeak::MovingPeak(ParamMap &v):Problem((v[param_proId]),(v[param_numDim]),string(),1),\
	DynamicContinuous((v[param_proId]), (v[param_numDim]),(v[param_numPeak]),1){
	setDimensionChange((v[param_flagNumDimChange]));
	setNumPeakChangeMode((v[param_peakNumChangeMode]));
	setNumPeaksChange((v[param_flagNumPeakChange]));
	setNoiseFlag((v[param_flagNoise]));
	setTimeLinkageFlag((v[param_flagTimeLinkage]));
	setChangeFre((v[param_changeFre]));
	setChangeType(CT_Random);
	
	m_vlength=(v[param_shiftLength]);
	if((v[param_flagNoise])==true)		setNoiseSeverity_((v[param_noiseSeverity]));
	if((v[param_flagTimeLinkage])==true)	setTimeLinkageSeverity((v[param_timelinkageSeverity]));
		
	if(!readData()) exit(0);
	m_peaksFound=0;
	allocateMemory(m_numDim,m_numPeaks);
	m_name="DYN_CONT_MovingPeak";
	initialize();
	setNumberofChanges((double)v[param_changeRatio]);
	m_globalOpt.setFlagLocTrue();
	m_proPar<<"Vlength:"<<m_vlength<<"; ";
	
}
MovingPeak::MovingPeak(const int rId, const int rDimNumber,  const int rNumPeaks,double const rChangingRatio,const bool rFlagDimChange,\
	const bool rFlagNumPeakChange,const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage):\
	Problem(rId,rDimNumber,string(),1),DynamicContinuous(rId,rDimNumber,rNumPeaks,1){

	setDimensionChange(rFlagDimChange);
	setNumPeakChangeMode(peakNumChangeMode);
	setNumPeaksChange(rFlagNumPeakChange);
    setNoiseFlag(flagNoise);
	setTimeLinkageFlag(flagTimelinkage);

	if(!readData()) exit(0);
	m_peaksFound=0;

	allocateMemory(m_numDim,m_numPeaks);
	m_name="DYN_CONT_MovingPeak";
	initialize();
	setNumberofChanges((int)(rChangingRatio*m_numPeaks));
	m_globalOpt.setFlagLocTrue();
	m_proPar<<"Vlength:"<<m_vlength<<"; ";


}
void MovingPeak::initialize(){
    int i=0;
	setAccuracy(0.1);
	setDisAccuracy(0.2);
	setSearchRange(0,100);
    setOptType(MAX_OPT);
	m_globalOpt.setRecordFlag(false);
	updateTimeLinkage();
	for ( i=0; i< m_numPeaks; i++)
		for ( int j=0;j<m_numDim; j++){
	  mpp_peak[i][j] = 100.0*Global::msp_global->mp_uniformPro->Next();
	  mpp_prevMovement[i][j] = Global::msp_global->mp_uniformPro->Next()-0.5;
	}

	if (m_standardHeight <= 0.0){
        for ( i=0; i< m_numPeaks; i++) mp_height[i]=(m_maxHeight-m_minHeight)*Global::msp_global->mp_uniformPro->Next()+m_minHeight;
	}else{
        for (i=0; i< m_numPeaks; i++) mp_height[i]= m_standardHeight;
	}

	if (m_standardWidth <= 0.0){
        for (i=0; i< m_numPeaks; i++)
            mp_width[i]= (m_maxWidth-m_minWidth)*Global::msp_global->mp_uniformPro->Next()+m_minWidth;
	}else{
	for (i=0; i< m_numPeaks; i++)
	  mp_width[i]= m_standardWidth;
	}

	calculateGlobalOptima();
	/*for (i=0; i< m_numPeaks; i++) {
		mp_heightOrder[i]=i;
		mp_found[i]=false;
	}
	vector<int> idx(m_numPeaks);
	gQuickSort(mp_height,m_numPeaks,idx);
	copy(idx.begin(),idx.end(),mp_heightOrder);
	gAmendSortedOrder<double*>(mp_height,mp_heightOrder,mp_amendedHeightOrder,m_numPeaks);*/
	for ( i=0; i< m_numPeaks; i++) mp_isTracked[i]=0;

	for (i=0;i<m_numPeaks; i++) 
	copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);
	//calculateAssociateRadius();
	m_peakQaulity=0;
	addProTag(MMP);
}
void MovingPeak::reset(){

	setSeverity();

    int i=0;
    for ( i=0; i< m_numPeaks; i++)
	for (int j=0;j<m_numDim; j++){
	  mpp_peak[i][j] = 100.0*Global::msp_global->mp_uniformPro->Next();
	  mpp_prevMovement[i][j] = Global::msp_global->mp_uniformPro->Next()-0.5;
	}

	if (m_standardHeight <= 0.0){
        for ( i=0; i< m_numPeaks; i++) mp_height[i]=(m_maxHeight-m_minHeight)*Global::msp_global->mp_uniformPro->Next()+m_minHeight;
	}else{
        for (i=0; i< m_numPeaks; i++) mp_height[i]= m_standardHeight;
	}

	if (m_standardWidth <= 0.0){
        for (i=0; i< m_numPeaks; i++)
            mp_width[i]= (m_maxWidth-m_minWidth)*Global::msp_global->mp_uniformPro->Next()+m_minWidth;
	}else{
	for (i=0; i< m_numPeaks; i++)
	  mp_width[i]= m_standardWidth;
	}

	calculateGlobalOptima();
    m_changeCounter=0;
	for ( i=0; i< m_numPeaks; i++) mp_found[i]=false;

	

	for (i=0;i<m_numPeaks; i++) 	copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	//calculateAssociateRadius();
	/*for (i=0; i< m_numPeaks; i++) {
			mp_heightOrder[i]=i;
			mp_found[i]=false;
	}
	vector<int> idx(m_numPeaks);
	gQuickSort(mp_height,m_numPeaks,idx);
	copy(idx.begin(),idx.end(),mp_heightOrder);

	
	gAmendSortedOrder<double*>(mp_height,mp_heightOrder,mp_amendedHeightOrder,m_numPeaks);*/
	m_peakQaulity=0;
	m_peaksFound=0;
}

/* free disc space at end of program */
MovingPeak::~MovingPeak()
{
    freeMemory();
    
}

/* current_peak_calc determines the peak of the current best individual */
void MovingPeak::currentPeakCalc(const double *gen)
{
  int i;
  double maximum = -100000.0, dummy;

  m_currentPeak = 0;
  maximum = functionSelection(gen, 0);
  for(i=1; i<m_numPeaks; i++){
    dummy = functionSelection(gen, i);
    if (dummy > maximum){
      maximum = dummy;
      m_currentPeak = i;
    }
  }
}

ReturnFlag MovingPeak::evaluate_(VirtualEncoding  &ss, bool rFlag, ProgramMode mode, bool flag2){
	CodeVReal &s=dynamic_cast<CodeVReal &>(ss);
	double *x=new double[m_numDim];
	copy(s.m_x.begin(),s.m_x.end(),x);
	if(this->m_noiseFlag)	addNoise(x);

    double maximum = LONG_MIN, dummy;

    for(int i=0; i<m_numPeaks; i++){
		//if(maximum>mp_height[i]) continue; //optimization on the obj evaluation
		dummy = functionSelection(x, i);
        if (dummy > maximum)      maximum = dummy;
    }

    if (m_useBasisFunction){
		dummy = functionSelection(x,-1);
		/* If value of basis function is higher return it */
		if (maximum < dummy)     maximum = dummy;
    }
	s.m_obj[0]=maximum;

	if(rFlag&&m_evals%m_changeFre==0){
		Solution<CodeVReal>::initilizeWB(s);
	}

	if(rFlag&&isTracked(x,s.m_obj)) updatePeakQaulity();
    if(rFlag)    m_evals++;
	bool flag;
	#ifdef OFEC_CONSOLE
		if(Global::msp_global->mp_algorithm!=nullptr)	flag=!Global::msp_global->mp_algorithm->ifTerminating();
		else flag=true;
	#endif
	#ifdef OFEC_DEMON
		flag=true;
	#endif
    if(rFlag&&m_evals%m_changeFre==0&&flag){
		
		//g_mutexStream.lock();
		//cout<<"The number of changes: "<<m_changeCounter<<endl;
		//g_mutexStream.unlock();
		//for(int i=0;i<m_numPeaks;i++)		printPeak(i);
		DynamicProblem::change();
		//for(int i=0;i<m_numPeaks;i++)		printPeak(i);
		//getchar();
	}
	delete [] x;
	x=0;
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


/* dummy evaluation function allows to evaluate without being counted */
double MovingPeak::dummyEval (const double *gen)
{
  int i;
  double maximum = -100000.0, dummy;

  for(i=0; i<m_numPeaks; i++){
    dummy = functionSelection(gen, i);
    if (dummy > maximum)    maximum = dummy;
    }

  if (m_useBasisFunction){
    dummy = functionSelection(gen,-1);
    /* If value of basis function is higher return it */
    if (maximum < dummy)    maximum = dummy;
  }
  return(maximum);
}

/* whenever this function is called, the peaks are changed */
void MovingPeak::randomChange(){
	int i=0,j=0;
	double sum, sum2, offset;

	for (i=0;i<m_numPeaks; i++) 	copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	for(i=0; i<m_numPeaks; i++){
		if(mp_whetherChange[i]==false) continue;
		/* shift peak locations */
		sum = 0.0;
		for ( j=0; j<m_numDim; j++){
			mp_shift[j]=Global::msp_global->mp_uniformPro->Next()-0.5;
			sum += mp_shift[j]*mp_shift[j];
		}
		if(sum>0.0)		sum = m_vlength/sqrt(sum);
		else  sum = 0.0;                        /* only in case of rounding errors */

		sum2=0.0;
		for (j=0; j<m_numDim; j++){
			mp_shift[j]=sum*(1.0-m_lambda)*mp_shift[j]+m_lambda*mpp_prevMovement[i][j];
			sum2 += mp_shift[j]*mp_shift[j];
		}
		if(sum2>0.0)sum2 = m_vlength/sqrt(sum2);
		else     sum2 = 0.0;                      /* only in case of rounding errors */

		for(j=0; j<m_numDim; j++){
			mp_shift[j]*=sum2;
			mpp_prevMovement[i][j]= mp_shift[j];
			if (m_searchRange[j].m_lower>(mpp_peak[i][j]+mpp_prevMovement[i][j]) ){
				mpp_peak[i][j] = 2.0*m_searchRange[j].m_lower-mpp_peak[i][j]-mpp_prevMovement[i][j];
				mpp_prevMovement[i][j]*=-1.0;
			}
			else if(m_searchRange[j].m_upper<(mpp_peak[i][j]+mpp_prevMovement[i][j]) ){
				mpp_peak[i][j]= 2.0*m_searchRange[j].m_upper-mpp_peak[i][j]-mpp_prevMovement[i][j];
				mpp_prevMovement[i][j]*=-1.0;
			}else
			mpp_peak[i][j] += mpp_prevMovement[i][j];
		}
		
		/* change peak width */

		offset = Global::msp_global->mp_normalPro->Next()*mp_widthSeverity[i];
		if ((mp_width[i]+offset) < m_minWidth)		mp_width[i] = 2.0*m_minWidth-mp_width[i]-offset;
		else if ((mp_width[i]+offset) > m_maxWidth)	mp_width[i]= 2.0*m_maxWidth-mp_width[i]-offset;
		else	mp_width[i] += offset;
		
		if(getChangeCounter()>1&&m_changePeakRatio<1.0&&isGOpt(i)) continue; 
		/* change peak height */

		offset = mp_heightSeverity[i]*Global::msp_global->mp_normalPro->Next();

		if ((mp_height[i]+offset) < m_minHeight)	mp_height[i] = 2.0*m_minHeight-mp_height[i]-offset;
		else if ((mp_height[i]+offset) > m_maxHeight)	mp_height[i]= 2.0*m_maxHeight-mp_height[i]-offset;
		else	mp_height[i] += offset;
	}

	calculateGlobalOptima();
	updateNumberofChanges();

}


/* Basis Functions */

/* This gives a constant value back to the eval-function that chooses the max of them */
double MovingPeak::constantBasisFunc(const double *gen)
{
  return 0.0;
}

double MovingPeak::fivePeakBasisFunc(const double *gen){
	double maximum = -100000.0, dummy=0;
	for(int i=0; i<5; i++){
		dummy = (gen[0]-basis_peak[i][0])*(gen[0]-basis_peak[i][0]);
		for (int j=1; j< m_numDim; j++)  dummy += (gen[j]-basis_peak[i][j])*(gen[j]-basis_peak[i][j]);
		dummy = basis_peak[i][m_numDim+1]-(basis_peak[i][m_numDim]*dummy);
		if (dummy > maximum)       maximum = dummy;
	}
	return maximum;
}

/* Peak Functions */

/* sharp peaks */
double MovingPeak::peakFunction1(const double *gen, int peak_number){
 
  double dummy = (gen[0]-mpp_peak[peak_number][0])*(gen[0]-mpp_peak[peak_number][0]);
  for (int j=1; j< m_numDim; j++)
    dummy += (gen[j]-mpp_peak[peak_number][j])*(gen[j]-mpp_peak[peak_number][j]);

  return mp_height[peak_number]/(1+mp_width[peak_number]*dummy);
}

double MovingPeak::peakFunctionHilly(const double *gen, int peak_number){
  int j=0;
  double dummy =  (gen[0]-mpp_peak[peak_number][0])*(gen[0]-mpp_peak[peak_number][0]);
  for (j=1; j< m_numDim; j++)
    dummy += (gen[j]-mpp_peak[peak_number][j])*(gen[j]-mpp_peak[peak_number][j]);

  return mp_height[peak_number]- mp_width[peak_number]*dummy-0.01*sin(20.0*dummy);
}

double MovingPeak::peakFunctionTwin(const double  *gen, int peak_number) /* two twin peaks moving together */
{
  int j;
  double maximum = -100000.0, dummy= pow(gen[0]-mpp_peak[peak_number][0],2);
  for (j=1; j< m_numDim; j++)
     dummy += pow(gen[j]-mpp_peak[peak_number][j],2);

  dummy = mp_height[peak_number]- mp_width[peak_number]*dummy;

  maximum = dummy;
  dummy = pow(gen[0]-(mpp_peak[peak_number][0]+twin_peak[0]),2);
  for (j=1; j< m_numDim; j++)
     dummy += pow(gen[j]-(mpp_peak[peak_number][j]+twin_peak[0]),2);

  dummy =  mp_height[peak_number]+twin_peak[m_numDim+1]-((mp_width[peak_number]+twin_peak[m_numDim])*dummy);
  if (dummy > maximum)
    maximum = dummy;

  return maximum;
}


double MovingPeak::peakFunctionCone(const double *gen, const int &peak_number){

	double val,  dummy=0;
	for (int j=0; j<m_numDim; j++){
		val=gen[j]-mpp_peak[peak_number][j];
		dummy += val*val;
	}
	if(dummy!=0)  dummy =mp_height[peak_number]-mp_width[peak_number]*sqrt(dummy);
	return dummy;
}
double MovingPeak::functionSelection(const double  *gen, const int &peak_number){
	double dummy=0;
	switch(m_F){
	case 1: {
			dummy=constantBasisFunc(gen);
			break;
		 }
	case 2: {
			dummy=fivePeakBasisFunc(gen);
			break;
		 }
		case 3: {
			dummy=peakFunction1(gen, peak_number);
			break;
		 }
	case 4: {
			dummy=peakFunctionCone(gen, peak_number);
			break;
		 }
	case 5: {
			dummy=peakFunctionHilly(gen, peak_number);
			break;
		 }
	case 6: {
			dummy=peakFunctionTwin(gen, peak_number);
			break;
		 }
	}
	return dummy;
}
/* The following procedures may be used to change the step size over time */


void MovingPeak::changeStepsizeRandom () /* assigns vlength a value from a normal distribution */
{
  m_vlength = Global::msp_global->mp_normalPro->Next();
}

void MovingPeak::changeStepsizeLinear() /* sinusoidal change of the stepsize, */
{
static	boost::thread_specific_ptr< int> counter;
if(!counter.get()) counter.reset(new int(1));
boost::thread_specific_ptr<double> frequency;
if(!frequency.get()) frequency.reset(new double (3.14159/20.0));

  m_vlength = 1+ sin((double)(*counter)*(*frequency));
  (*counter) ++;
}

int MovingPeak::getRightPeak()  /* returns 1 if current best individual is on highest mpp_peak, 0 otherwise */
{
    bool flag=false;

    for(int i=0;i<m_numPeaks;i++){
        if(mp_globalOptimaIdx[i]==true && m_currentPeak == i){
            flag=true;
            break;
        }
    }

    return flag;
}
void MovingPeak::setVlength(const double s){
	m_vlength=s;

    size_t start, end=0;
    start=m_proPar.str().find("Vlength:");
    for(auto i=start;i<m_proPar.str().size();i++){
        if(m_proPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Vlength:"<<m_vlength<<"; ";
    string result=m_proPar.str();
    result.replace(start,end-start+1, ss.str());
    m_proPar.str(result);

}

void MovingPeak::changeNumPeaks(){

	MovingPeak* mpb=new MovingPeak(Global::msm_pro["DYN_CONT_MovingPeak"],m_numDim,m_numPeaksTemp,m_changePeakRatio,m_flagDimensionChange
		,m_flagNumPeaksChange,m_numPeaksChangeMode,m_noiseFlag,m_timeLinkageFlag);

	mpb->parameterSetting(this);
	mpb->calculateGlobalOptima();

	freeMemory();
	DynamicContinuous::freeMemory();

	allocateMemory(m_numDim,m_numPeaksTemp);
    DynamicContinuous::allocateMemory(m_numDim,m_numPeaksTemp);
    m_numPeaks=m_numPeaksTemp;
	*this=*mpb;
	delete mpb;
	mpb=0;
}

void MovingPeak::parameterSetting(Problem * rP){
	DynamicContinuous::parameterSetting(rP);

	MovingPeak *mpb=dynamic_cast<MovingPeak *>(rP);
	unsigned dim=m_dimNumberTemp<mpb->getNumDim()?m_dimNumberTemp:mpb->getNumDim();
	int peaks=m_numPeaks<mpb->getNumberofPeak()?m_numPeaks:mpb->getNumberofPeak();

	m_F=mpb->m_F;
	m_vlength=mpb->m_vlength;
	m_lambda=mpb->m_lambda;
	m_useBasisFunction=mpb->m_useBasisFunction;
	m_calculateRightPeak=mpb->m_calculateRightPeak ;
	m_standardHeight=mpb->m_standardHeight;
	m_standardWidth=mpb->m_standardWidth;

	//gCopy(mp_shift,mpb->mp_shift,dim);
	copy(mpb->mp_shift,mpb->mp_shift+dim,mp_shift);
	//gCopy(mp_coveredPeaks,mpb->mp_coveredPeaks,peaks);
	copy(mpb->mp_coveredPeaks,mpb->mp_coveredPeaks+peaks,mp_coveredPeaks);

	for (int i=0; i<peaks; i++){
		//gCopy(mpp_prevMovement[i],mpb->mpp_prevMovement[i],dim);
		copy(mpb->mpp_prevMovement[i],mpb->mpp_prevMovement[i]+dim,mpp_prevMovement[i]);
	}

}
MovingPeak &MovingPeak::operator=(MovingPeak &other){
	if(this==&other) return *this;

	if(m_numDim!=other.m_numDim||m_numPeaks!=other.m_numPeaks) throw myException("Moving Peak assignment@MovingPeak::operator=");
	DynamicContinuous::operator=(other);

	m_F=other.m_F;
	m_vlength=other.m_vlength;
	m_lambda=other.m_lambda;
	m_useBasisFunction=other.m_useBasisFunction;
	m_calculateRightPeak=other.m_calculateRightPeak ;
	m_standardHeight=other.m_standardHeight;
	m_standardWidth=other.m_standardWidth;

	copy(other.mp_shift,other.mp_shift+m_numDim,mp_shift);
	copy(other.mp_coveredPeaks,other.mp_coveredPeaks+m_numPeaks,mp_coveredPeaks);


	for (int i=0; i<m_numPeaks; i++){
		copy(other.mpp_prevMovement[i],other.mpp_prevMovement[i]+m_numDim,mpp_prevMovement[i]);

	}
	return *this;
}

void MovingPeak::updatePeakQaulity(){
	/*relative height over the heightest; Note perform just when a peak is found */
	m_peakQaulity=0;
	double sum=0;
	for(int i=0;i<m_numPeaks;i++){
		if(mp_found[i]) m_peakQaulity+=mp_amendedHeightOrder[i]*1./mp_amendedHeightOrder[mp_heightOrder[m_numPeaks-1]];
		sum+=mp_amendedHeightOrder[i]*1./mp_amendedHeightOrder[mp_heightOrder[m_numPeaks-1]];
 	}
	if(sum>0)	m_peakQaulity/=sum;
	else m_peakQaulity=0;
}


void MovingPeak::calculateAssociateRadius(){
	// to calculate an assosiate radius of peak i, find the nearest peak j, get the valley point, the distance
	// between peak i and the valley point is the associate radius of peak i
	double *point=new double[m_numDim];
	double *asRa=new double[m_numDim];

	for(int i=0;i<m_numPeaks;i++){
		double dis; int nearest=-1;
		if(!isVisable(i)) continue;
		for(int j=0,count=0;j<m_numPeaks;j++,count++){
			if(j==i||!isVisable(j)) {count--;continue;}
			double d=0;
			for(int dim=0;dim<m_numDim;dim++){
				d+=(mpp_peak[i][dim]-mpp_peak[j][dim])*(mpp_peak[i][dim]-mpp_peak[j][dim]);
			}
			d=sqrt(d);
			if(0==count){
				dis=d;
				nearest=j;
			}else if(dis>d){
				dis=d;
				nearest=j;
			}
		}
		if(nearest!=-1){
			//normalize vector point
			for(int dim=0;dim<m_numDim;dim++){
				point[dim]=mpp_peak[nearest][dim]-mpp_peak[i][dim];
				point[dim]/=dis;
			}


			double height,asHeight;
			height=asHeight=mp_height[i];
			copy(mpp_peak[i],mpp_peak[i]+m_numDim,asRa);
			//test in direction of point with a step of dis/100
			while(asHeight<=height){
				bool flagBreak=false;
				for(int dim=0;dim<m_numDim;dim++){
					asRa[dim]+=dis/100*point[dim];
					if((asRa[dim]-mpp_peak[i][dim])*(mpp_peak[nearest][dim]-asRa[dim])<0){
						flagBreak=true;
						break;
					}
				}
				if(flagBreak) break;
				height=asHeight;
				CodeVReal s(m_numDim,m_numObj);
				copy(asRa,asRa+m_numDim,s.m_x.begin());
				evaluate_(s,false);
				asHeight=s.m_obj[0];
			}

			mp_associateRadius[i]=0;
			for(int dim=0;dim<m_numDim;dim++){
				// correction for one step backward
				asRa[dim]-=dis/200*point[dim];
				mp_associateRadius[i]+=(asRa[dim]-mpp_peak[i][dim])*(asRa[dim]-mpp_peak[i][dim]);
			}
			mp_associateRadius[i]=sqrt(mp_associateRadius[i]);
		}else{
			double *r=new double[2*m_numDim];
			for(int dim=0;dim<m_numDim;dim++){
				double u,l;
				m_searchRange.getSearchRange(l,u,dim);
				r[dim*2]=fabs(l-mpp_peak[i][dim]);
				r[dim*2+1]=fabs(u-mpp_peak[i][dim]);
			}
			mp_associateRadius[i]=*min_element(r,r+2*m_numDim);
			delete []r;
			r=0;
		}

	}

		delete [] point;
		delete [] asRa;
		point=0;asRa=0;
}
double MovingPeak::getVLength(){
	return m_vlength;
}
