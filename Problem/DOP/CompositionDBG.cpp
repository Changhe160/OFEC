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

#include "CompositionDBG.h"
#include "../../Global/global.h"
boost::thread_specific_ptr<ComDBGFuncID> CompositionDBG::ms_funID;

void CompositionDBG::freeMemory(){
    delete [] mp_stretchSeverity; 
	mp_stretchSeverity=0;
	delete [] mp_convergeSeverity;
	mp_convergeSeverity=0;

	delete [] mp_comFunction;
	mp_comFunction=0;
	delete []mp_fmax;
	mp_fmax=0;
}

CompositionDBG::CompositionDBG(ParamMap &v):Problem((v[param_proId]),(v[param_numDim]),string(),1),\
	RealDBG(v[param_proId],v[param_numDim],v[param_numPeak]),mv_comBoundary(msc_numComFuns){

    allocateMemory(m_numDim,m_numPeaks);
    m_name="DYN_CONT_CompositionDBG";
    setComBoundary();
	m_heightNormalizeSeverity=2000.;

	initialize(static_cast<ChangeType>((int)(v[param_changeType])),static_cast<ComDBGFuncID>((int)(v[param_comDBGFunID])),\
		(v[param_changeRatio]),(v[param_flagNumDimChange]),(v[param_flagNumPeakChange]),\
		(v[param_peakNumChangeMode]),(v[param_flagNoise]),(v[param_flagTimeLinkage]));
}

CompositionDBG::CompositionDBG(const int rId, const int rDimNumber, const int rNumPeaks,const ChangeType rT, \
							   const ComDBGFuncID rF, const double rChangingRatio,const bool rFlagDimChange,const bool rFlagNumPeakChange,\
							   const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage):\
							   Problem(rId,rDimNumber,string(),1),RealDBG(rId,rDimNumber,rNumPeaks),mv_comBoundary(msc_numComFuns){

    allocateMemory(m_numDim,m_numPeaks);
    m_name="DYN_CONT_CompositionDBG";
    setComBoundary();
	m_heightNormalizeSeverity=2000.;
	
	initialize(rT,rF,rChangingRatio,rFlagDimChange,rFlagNumPeakChange, peakNumChangeMode, flagNoise,flagTimelinkage);

}

void CompositionDBG::allocateMemory(const int rDimNum, const int rPeaks){

	
	for(int i=0;i<msc_numComFuns;i++){ 
		mv_comBoundary[i].resizeDim(rDimNum);
	}

	mp_convergeSeverity =new double[rPeaks];
	mp_stretchSeverity= new double[rPeaks];
	mp_comFunction=new ProblemTag[rPeaks];
	mp_fmax= new double[rPeaks];

}

CompositionDBG::~CompositionDBG(){
    freeMemory();
}
void CompositionDBG::setComBoundary(){
	//Sphere=0,Rastrigin,Weierstrass,Griewank,Ackley

	for(int j=0;j<m_numDim;j++){
			mv_comBoundary[0][j].m_upper=100.;
			mv_comBoundary[0][j].m_lower=-100.;
			mv_comBoundary[1][j].m_upper=5.;
			mv_comBoundary[1][j].m_lower=-5.;
			mv_comBoundary[2][j].m_upper=0.5;
			mv_comBoundary[2][j].m_lower=-0.5;
			mv_comBoundary[3][j].m_upper=100.;
			mv_comBoundary[3][j].m_lower=-100.;
			mv_comBoundary[4][j].m_upper=32.;
			mv_comBoundary[4][j].m_lower=-32.;
		}
}
void CompositionDBG::getComBoundary(const ProblemTag &f,double &l,double &u,const int rDimIdx)const{
	switch(f){
		case Sphere:
			l=mv_comBoundary[0][rDimIdx].m_lower;
			u=mv_comBoundary[0][rDimIdx].m_upper;
			break;
		case Rastrigin:
			l=mv_comBoundary[1][rDimIdx].m_lower;
			u=mv_comBoundary[1][rDimIdx].m_upper;
			break;
		case Weierstrass:
			l=mv_comBoundary[2][rDimIdx].m_lower;
			u=mv_comBoundary[2][rDimIdx].m_upper;
			break;
		case Griewank:
			l=mv_comBoundary[3][rDimIdx].m_lower;
			u=mv_comBoundary[3][rDimIdx].m_upper;
			break;
		case Ackley:
			l=mv_comBoundary[4][rDimIdx].m_lower;
			u=mv_comBoundary[4][rDimIdx].m_upper;
			break;
		default:
			throw myException("No the function in the basic component fs@CompositionDBG::getComBoundary");
			break;
	}
}
CompositionDBG & CompositionDBG::operator=(const CompositionDBG &rCom){
	if(this==&rCom) return *this;
	if(m_numDim!=rCom.m_numDim||m_numPeaks!=rCom.m_numPeaks) throw myException("CompositionDBG assignment"); 
	RealDBG::operator =(rCom);

	setComBoundary();
	copy(rCom.mp_convergeSeverity,rCom.mp_convergeSeverity+m_numPeaks,mp_convergeSeverity);
	copy(rCom.mp_stretchSeverity,rCom.mp_stretchSeverity+m_numPeaks,mp_stretchSeverity);

	m_heightNormalizeSeverity=rCom.m_heightNormalizeSeverity;
	copy(rCom.mp_comFunction,rCom.mp_comFunction+m_numPeaks,mp_comFunction);
    copy(rCom.mp_fmax,rCom.mp_fmax+m_numPeaks,mp_fmax);

	return *this;
}
void CompositionDBG::setRotationMatrix(){
	// for each basic function of dimension n(even number), R=R(l1,l2)*R(l3,l4)*....*R(ln-1,ln), 0<=li<=n
	Matrix I(m_numDim,m_numDim);

	int * d=new int[m_numDim];
	Global::msp_global->initializeRandomArray(d,m_numDim,Program_Problem);
	for(int i=0;i<m_numPeaks;i++){
		for(int j=0;j+1<m_numDim;j+=2){
			double angle=2*OFEC_PI*Global::msp_global->mp_uniformPro->Next();				// random angle for rotation plane of d[j]-d[j+1] from d[j]th axis to d[j+1]th axis
			I.setRotationAngle(d[j],d[j+1],angle);
			if(j==0)  mp_rotationMatrix[i]=I;
			else	mp_rotationMatrix[i]*=I;
			I.identity();
		}
	}
	delete [] d;
	d=0;
}
void CompositionDBG::setCovergeSevrity( const double *cs){
	copy(cs,cs+m_numPeaks,mp_convergeSeverity);
}
void CompositionDBG::setStretchSeverity(){

	for(int i=0;i<m_numPeaks;i++){
		double l,u;
        getComBoundary(mp_comFunction[i],l,u);
        mp_stretchSeverity[i]=mp_convergeSeverity[i]*(m_searchRange[0].m_upper-m_searchRange[0].m_lower)/(u-l);
	}
}
void CompositionDBG::setBasicFunction(const ProblemTag *bf){
	copy(bf,bf+m_numPeaks,mp_comFunction);
}
ReturnFlag CompositionDBG::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag_){
	CodeVReal &s=dynamic_cast<CodeVReal &>(ss);
	double *x=new double[m_numDim];
	copy(s.m_x.begin(),s.m_x.end(),x);

	if(this->m_noiseFlag)	addNoise(x);
	
	vector<double> width(m_numPeaks,0),fit(m_numPeaks);
	for(int i=0;i<m_numPeaks;i++){ // calculate weight for each function		
		for(int j=0;j<m_numDim;j++)
			width[i]+=(x[j]-mpp_peak[i][j])*(x[j]-mpp_peak[i][j]);
		if(width[i]!=0)	width[i]=exp(-sqrt(width[i]/(2*m_numDim*mp_convergeSeverity[i]*mp_convergeSeverity[i])));
	}

	for(int i=0;i<m_numPeaks;i++){ // calculate objective value for each function
		for(int j=0;j<m_numDim;j++)	// calculate the objective value of tranformation function i
			x[j]=(x[j]-mpp_peak[i][j])/mp_stretchSeverity[i];//((1+fabs(mpp_peak[i][j]/mp_searchRange[j].m_upper))*
		Matrix m(m_numDim,1);
		m.setDataRow(x,m_numDim);
		m*=mp_rotationMatrix[i];
		copy(m[0].begin(),m[0].end(),x);
		correctSolution(mp_comFunction[i],x);
		fit[i]=selectFun(mp_comFunction[i],x);
		fit[i]=m_heightNormalizeSeverity*fit[i]/fabs(mp_fmax[i]);
		copy(s.m_x.begin(),s.m_x.end(),x);
	}
	double sumw=0,wmax;
	wmax=*max_element(width.begin(),width.end());
	for(int i=0;i<m_numPeaks;i++)
		if(width[i]!=wmax)
			width[i]=width[i]*(1-pow(wmax,10));
	for(int i=0;i<m_numPeaks;i++)
		sumw+=width[i];
	for(int i=0;i<m_numPeaks;i++)
		width[i]/=sumw;
	double obj=0;
	for(int i=0;i<m_numPeaks;i++)
		obj+=width[i]*(fit[i]+mp_height[i]);
	s.m_obj[0]=obj;

	if(rFlag&&m_evals%m_changeFre==0) Solution<CodeVReal>::initilizeWB(s);

	if(rFlag){
		isTracked(x,s.m_obj);
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
    if(rFlag&&m_evals%m_changeFre==0&&flag) {
		DynamicProblem::change();
		if(m_timeLinkageFlag) updateTimeLinkage();
	}
	delete []x;
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

inline double CompositionDBG::selectFun(const ProblemTag &f,double *x){
	double value;
	switch(f){
	case Sphere:
		value=fSphere(x);
		break;
	case Rastrigin:
		value=fRastrigin(x);
		break;
	case Weierstrass:
		value=fWeierstrass(x);
		break;
	case Griewank:
		value=fGriewank(x);
		break;
	case Ackley:
		value=fAckley(x);
		break;
	default:
		break;
	}
	return value;
}
inline double CompositionDBG::fAckley(double *x){
	double fitness=0;
	double s1=0,s2=0;
	for(int i=0;i<m_numDim;i++){
		s1+=x[i]*x[i];
		s2+=cos(2*OFEC_PI*x[i]);
	}
	fitness=-20*exp(-0.2*sqrt(s1/m_numDim))-exp(s2/m_numDim)+20+OFEC_E;
	return fitness;
}
inline double CompositionDBG::fGriewank(double *x){
	double s1=0,s2=1;
	for(int i=0;i<m_numDim;i++){
		s1+=x[i]*x[i]/4000.;
		s2*=cos(x[i]/sqrt((double) (i+1)));
	}
	return s1-s2+1.;
}
inline double CompositionDBG::fRastrigin(double *x){
	double fit=0;
	for(int i=0;i<m_numDim;i++)
		fit=fit+x[i]*x[i]-10.*cos(2*OFEC_PI*x[i])+10.;
	return fit;
}
inline double CompositionDBG::fSphere(double *x){
	double fit=0;
	for(int i=0;i<m_numDim;i++)
		fit+=x[i]*x[i];
	return fit;
}
inline double CompositionDBG::fWeierstrass(double *x){
	double a=0.5,b=3;
	int kmax=20;
	double fit=0,s=0;
	for(int i=0;i<m_numDim;i++)
		for(int k=0;k<=kmax;k++)
			fit+=pow(a,k)*cos(2*OFEC_PI*pow(b,k)*(x[i]+0.5));
	for(int k=0;k<=kmax;k++)
			s+=pow(a,k)*cos(2*OFEC_PI*pow(b,k)*0.5);
	s=s*m_numDim;
	return fit-s;
}
void CompositionDBG::correctSolution(const ProblemTag &f, double *x){

	double l,u;
	getComBoundary(f,l,u);
	for(int j=0;j<m_numDim;j++){
		if(x[j]>u)  x[j]=u;
		else if(x[j]<l)  x[j]=l;
	}
}

void CompositionDBG::randomChange(){
    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	//change the global minimum value of each function
	heightStandardChange();
	//change the position of global optimum of each function randomly
	positionStandardChange(0);

    restoreInfor();
	calculateGlobalOptima();
    updateNumberofChanges();
	m_changeType.counter++;
}
void CompositionDBG::recurrentNoisyChange(){ 
	for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	double initial_angle;
	double height_range=m_maxHeight-m_minHeight;

	double noisy;
	for(int i=0;i<m_numPeaks;i++){
	    if(mp_whetherChange[i]==false) continue;
		initial_angle=(double)getPeriod()*i/m_numPeaks;

		 mp_height[i]=sinValueNoisy(m_changeType.counter,m_minHeight,m_maxHeight,height_range,initial_angle,m_noisySeverity);
	}

	initial_angle=OFEC_PI*(sin(2*OFEC_PI*(m_changeType.counter)/m_period)+1)/12.;
	noisy=m_noisySeverity*Global::msp_global->mp_normalPro->Next();
	positionStandardChange(initial_angle+noisy);

    restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void CompositionDBG::recurrentChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);
	double initial_angle;
	double height_range=m_maxHeight-m_minHeight;

	for(int i=0;i<m_numPeaks;i++){
	    if(mp_whetherChange[i]==false) continue;
		initial_angle=(double)m_period*i/m_numPeaks;
		mp_height[i]=m_minHeight+height_range*(sin(2*OFEC_PI*(m_changeType.counter+initial_angle)/m_period)+1)/2.;
	}
	initial_angle=OFEC_PI*(sin(2*OFEC_PI*m_changeType.counter/m_period)+1)/12.;
	positionStandardChange(initial_angle);

    restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void CompositionDBG::smallStepChange(){
    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	heightStandardChange();
	positionStandardChange(0);

    restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void CompositionDBG::largeStepChange(){
    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	heightStandardChange();
	positionStandardChange(0);

	restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}
void CompositionDBG::chaoticChange(){

    for (int i=0;i<m_numPeaks; i++) copy(mpp_peak[i],mpp_peak[i]+m_numDim,mpp_prePeak[i]);
	copy(mp_height,mp_height+m_numPeaks,mp_preHeight);
	copy(mp_width,mp_width+m_numPeaks,mp_preWidth);

	for(int i=0;i<m_numPeaks;i++){
	    if(mp_whetherChange[i]==false) continue;
		mp_height[i]=gChaoticValue(mp_height[i],m_minHeight,m_maxHeight);
	}

	positionStandardChange(0);

    restoreInfor();
	calculateGlobalOptima();
	 updateNumberofChanges();
	m_changeType.counter++;
}

void CompositionDBG::changeDimension(){
    /// no need to preserve  previous information, e.g., positions, heights, width....

	CompositionDBG* r_dbg=new CompositionDBG(m_id,m_dimNumberTemp,m_numPeaks,m_changeType.type,*ms_funID,m_changePeakRatio,m_flagDimensionChange,
		m_flagNumPeaksChange,m_numPeaksChangeMode,m_noiseFlag,m_timeLinkageFlag);


	if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy){
        r_dbg->setPeriod(m_period);
	}

	r_dbg->parameterSetting(this);
	r_dbg->calculateGlobalOptima();

	freeMemory();
	RealDBG::freeMemory();
	DynamicContinuous::freeMemory();
	Problem::freeMemory();

	Problem::allocateMemory(m_dimNumberTemp);
	ContinuousProblem::allocateMemory(m_dimNumberTemp);
    DynamicContinuous::allocateMemory(m_dimNumberTemp,m_numPeaks);
    RealDBG::allocateMemory(m_dimNumberTemp,m_numPeaks);
    allocateMemory(m_dimNumberTemp,m_numPeaks);

    m_numDim=m_dimNumberTemp;
    setPeriod(m_period);
	*this=*r_dbg;
	delete r_dbg;
	r_dbg=0;

}

void CompositionDBG::parameterSetting(Problem * rDP){
	RealDBG::parameterSetting(rDP);
	RealDBG *rP=dynamic_cast<RealDBG*>(rDP);

	CompositionDBG* r_dbg=static_cast<CompositionDBG*>(rP);

    int peaks=m_numPeaks<rP->getNumberofPeak()?m_numPeaks:rP->getNumberofPeak();

	copy(r_dbg->mp_comFunction,r_dbg->mp_comFunction+peaks,mp_comFunction);
//	setRotationMatrix();
	copy(r_dbg->mp_convergeSeverity,r_dbg->mp_convergeSeverity+peaks,mp_convergeSeverity);
	copy(r_dbg->mp_stretchSeverity,r_dbg->mp_stretchSeverity+peaks,mp_stretchSeverity);

}


void CompositionDBG::initialize(const ChangeType rT, const ComDBGFuncID rF,double const rChangingRatio,const bool rFlagDimChange, const bool rFlagNumPeakChange,
	const int peakNumChangeMode,const bool flagNoise, const bool flagTimelinkage){

    RealDBG::initialize(rT,rFlagDimChange,rFlagNumPeakChange,peakNumChangeMode,flagNoise,flagTimelinkage);
    setOptType(MIN_OPT);
    setNumberofChanges(rChangingRatio);
	if(!ms_funID.get()) ms_funID.reset(new ComDBGFuncID());
    *CompositionDBG::ms_funID=rF;
	setAccuracy(1.0);
	setDisAccuracy(0.1);
	m_globalOpt.setFlagLocTrue();
    ProblemTag *basic_fun=new ProblemTag[m_numPeaks];

	switch(*CompositionDBG::ms_funID){
		case COMDBG_SPHERE:
			for(int i=0;i<m_numPeaks;i++) basic_fun[i]=Sphere;
			break;
		case COMDBG_RASTRIGIN:
			for(int i=0;i<m_numPeaks;i++) basic_fun[i]=Rastrigin;
			break;
		case COMDBG_GRIEWANK:
			for(int i=0;i<m_numPeaks;i++) basic_fun[i]=Griewank;
			break;
		case COMDBG_ACKLEY:
			for(int i=0;i<m_numPeaks;i++) basic_fun[i]=Ackley;
			break;
		case COMDBG_HYBRID:
			basic_fun[0]=Sphere;		basic_fun[5]=Sphere;
			basic_fun[1]=Rastrigin;		basic_fun[6]=Rastrigin;
			basic_fun[2]=Weierstrass;	basic_fun[7]=Weierstrass;
			basic_fun[3]=Griewank;		basic_fun[8]=Griewank;
			basic_fun[4]=Ackley;		basic_fun[9]=Ackley;
			for(int i=10;i<m_numPeaks;i++) basic_fun[i]=Sphere;
			break;
	}
	setBasicFunction(basic_fun);
    delete []basic_fun;
	basic_fun=0;

	double *t=new double[m_numPeaks];
	for(int i=0;i<m_numPeaks;i++)t[i]=1.;
	setCovergeSevrity(t);
	setStretchSeverity();
	setRotationMatrix();
	
	Matrix m(m_numDim,1);
	double *gene=new double[m_numDim];
	for(int i=0;i<m_numPeaks;i++){
		for(int j=0;j<m_numDim;j++){ // calculate the estimate max value of funciton i
			gene[j]=m_searchRange[j].m_upper;
			gene[j]/=mp_stretchSeverity[i];
		}
		m.setDataRow(gene,m_numDim);
		m*=mp_rotationMatrix[i];
		copy(m[0].begin(),m[0].end(),gene);
		correctSolution(mp_comFunction[i],gene);
		mp_fmax[i]=selectFun(mp_comFunction[i],gene);
		if(mp_fmax[i]==0)   throw myException("the estimation max value must be greater not equal to 0@CompositionDBG::initialize");

	}
	
	calculateGlobalOptima();
	updateTimeLinkage();

	delete [] t;
	delete [] gene;
	t=0;
	gene=0;

 }

void  CompositionDBG::changeNumPeaks(){

    CompositionDBG* r_dbg=new CompositionDBG(m_id,m_numDim,m_numPeaksTemp,m_changeType.type,*ms_funID,
                                             m_changePeakRatio,m_flagDimensionChange,m_flagNumPeaksChange,
											 m_numPeaksChangeMode,m_noiseFlag,m_timeLinkageFlag);

	if(m_changeType.type==CT_Recurrent||m_changeType.type==CT_RecurrentNoisy)
		r_dbg->setPeriod(m_period);
	r_dbg->parameterSetting(this);

	r_dbg->calculateGlobalOptima();

	freeMemory();
	RealDBG::freeMemory();
	DynamicContinuous::freeMemory();

    DynamicContinuous::allocateMemory(m_numDim,m_numPeaksTemp);
    RealDBG::allocateMemory(m_numDim,m_numPeaksTemp);
    allocateMemory(m_numDim,m_numPeaksTemp);

    m_numPeaks=m_numPeaksTemp;
    setPeriod(m_period);
	*this=*r_dbg;
	delete r_dbg;
	r_dbg=0;
}



