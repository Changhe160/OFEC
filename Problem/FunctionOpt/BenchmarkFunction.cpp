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
// Created: 20 July 2011
// Last modified:
#include "BenchmarkFunction.h"

BenchmarkFunction::~BenchmarkFunction(){
    freeMemory();
}
BenchmarkFunction::BenchmarkFunction(const int rId, const int rDimNumber,   string rName, int numObj):\
	ContinuousProblem(rId,rDimNumber,rName,numObj), \
	m_scaleFlag(false),m_rotationFlag(false),m_translationFlag(false),m_biasFlag(false),m_scale(1.0),m_bias(0.),\
	m_originalGlobalOpt(rDimNumber,numObj),m_idx(0),m_gidx(0){
	
	setAccuracy(1.0e-6);
    allocateMemory(m_numDim);
    zeroTranslation();
    mp_rotationMatrix->identity();
}

void BenchmarkFunction::allocateMemory(const int rDimNumber ){
    mp_translation=new double[rDimNumber];
	mp_rotationMatrix =new Matrix(rDimNumber,rDimNumber);
	
}

void BenchmarkFunction::freeMemory(){
    if(mp_translation) delete [] mp_translation;

    if(mp_rotationMatrix)       delete mp_rotationMatrix;

    mp_translation=0;
    mp_rotationMatrix=0;
    

}

BenchmarkFunction& BenchmarkFunction::operator=(const BenchmarkFunction & rBF){
    if(this==&rBF) return *this;

    if(m_numDim!=rBF.m_numDim){
         throw myException("The number of dimensions must be same!@BenchmarkFunction::operator=");
    }

	ContinuousProblem::operator=(rBF);

    m_scaleFlag=rBF.m_scaleFlag;
    m_biasFlag=rBF.m_biasFlag;
    m_rotationFlag=rBF.m_rotationFlag;
    m_translationFlag=rBF.m_translationFlag;

    m_scale=rBF.m_scale;
    m_bias=rBF.m_bias;
    

	copy(rBF.mp_translation,rBF.mp_translation+m_numDim,mp_translation);
    *mp_rotationMatrix=*rBF.mp_rotationMatrix;

    m_conditionNumber=rBF.m_conditionNumber;

	//TODO debug
	m_originalGlobalOpt=rBF.m_originalGlobalOpt;
	m_idx = rBF.m_idx;
	m_gidx = rBF.m_gidx;
    return *this;
}

void BenchmarkFunction::setTranslation(const double *rTrans){
	copy(rTrans,rTrans+m_numDim,mp_translation);
    if(isZeroTranslation()) m_translationFlag=false;
    else m_translationFlag=true;

}
void BenchmarkFunction::setRotation(const double * const * rRot){
    mp_rotationMatrix->setData(rRot);

    if(mp_rotationMatrix->isIdentity()) m_rotationFlag=false;
    else m_rotationFlag=true;

}
void BenchmarkFunction::setBias(double rBias){
    m_bias=rBias;
    if(m_bias!=0) m_biasFlag=true;
    else m_biasFlag=false;
}
void BenchmarkFunction::setScale(double rScale){
    m_scale=rScale;

    if(m_scale!=1) m_scaleFlag=true;
    else m_scaleFlag=false;
}

void BenchmarkFunction::zeroTranslation(){
     for(int i=0;i<m_numDim;i++) mp_translation[i]=0.;

}

bool BenchmarkFunction::isZeroTranslation(){
    for(int i=0;i<m_numDim;i++){
        if(mp_translation[i]!=0) return false;
    }

    return true;

}

void BenchmarkFunction::transform(double * x){

    if(!m_translationFlag&&!m_scaleFlag&&!m_rotationFlag) return;


    double *x_=new double[m_numDim];
    // gCopy(x_,x,m_numDim);
	 copy(x,x+m_numDim,x_);
    if(m_translationFlag){
        for( int i=0;i<m_numDim; i++ ) 	{
             x_[ i ] = x[ i ] - mp_translation[ i ];
        }
	}

    if(m_scaleFlag){
        for( int i=0;i<m_numDim; i++ )     x_[ i ]/=m_scale;
    }

    if(m_rotationFlag){
        for(int i=0;i<m_numDim; i++ ) {
            x[i] = 0;

		for(int j = 0; j < m_numDim; j++ ) {
			x[ i ] += (*mp_rotationMatrix)[ j ][ i ] * x_[ j ];
		}
        }
    }
     else //gCopy(x,x_,m_numDim);
		 copy(x_,x_+m_numDim,x);

    delete [] x_;
	x_=0;
}
void BenchmarkFunction::setConditionNumber(double rC){
    m_conditionNumber=rC;
}

bool BenchmarkFunction::loadTranslation(){
    	// Initial the location of shifted global optimum
	if(!m_originalGlobalOpt.flagLoc()){
            throw myException("error, the original global optimia is not defined!@BenchmarkFunction::loadTranslation");
    }

	string s;
	stringstream ss;
	ss<<m_numDim<<"Dim.txt";
	s=ss.str();
	s.insert(0, m_name+"_Opt_");
	s.insert(0,"Problem/FunctionOpt/Data/" );//probDataPath
	s.insert(0,Global::g_arg[param_workingDir]);
	ifstream in;
	in.open(s.data());
	if(in.fail()){
		if(IS_PROBLEM_NAME(m_id,"FUN_RS_Ackley_Bound_CEC05")){
			for(int j=0;j<m_numDim;j++){
				if((j+1)%2==1)
					mp_translation[j]=m_searchRange[j].m_lower;
				else
					mp_translation[j]=m_searchRange[j].m_lower+(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*(1-Global::msp_global->mp_uniformPro->Next());
			}
		}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Schwefel")||IS_PROBLEM_NAME(m_id,"FUN_S_Schwefel")){

			for(int j=0;j<m_numDim;j++){

				mp_translation[j]=(Global::msp_global->mp_uniformPro->Next()-0.5)*20;
			}
		}
		else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank_noBounds_CEC05")){ 
			for(int j=0;j<m_numDim;j++){
				mp_translation[j] = -600.0 + Global::msp_global->mp_uniformPro->Next() * 600;
			}
		}else{
			for(int j=0;j<m_numDim;j++){
				double rl,ru,range,x;
				x=m_originalGlobalOpt[0].data().m_x[j];
				ru=m_searchRange[j].m_upper-x; 
				rl=x-m_searchRange[j].m_lower;
				range=rl<ru?rl:ru;
				mp_translation[j]=(Global::msp_global->mp_uniformPro->Next()-0.5)*2*range;
			}

		}

		ofstream out(s.c_str());
        for(int j=0;j<m_numDim;j++)        out<<mp_translation[j]<<" ";
		out.close();

	}else{

        for(int j=0;j<m_numDim;j++) {
            in>>mp_translation[j];
		}
	}
	in.close();
	m_translationFlag=true;
    return true;

}
bool BenchmarkFunction::loadRotation(){
	string s;
	stringstream ss;
	ss<<m_numDim<<"Dim.txt";
	s=ss.str();
	s.insert(0,m_name+"_RotM_");

	s.insert(0,"Problem/FunctionOpt/Data/" );//probDataPath
	s.insert(0,Global::g_arg[param_workingDir]);
	ifstream in;
	in.open(s.data());
	if(in.fail()){
		if(IS_PROBLEM_NAME(m_id,"FUN_RS_Elliptic_CEC05")){
            mp_rotationMatrix->identity();
		}else if(IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank")||IS_PROBLEM_NAME(m_id,"FUN_RS_Griewank_noBounds_CEC05")){
			mp_rotationMatrix->randomize(Program_Problem);
            mp_rotationMatrix->generateRotationMatrix(m_conditionNumber,Program_Problem);
			(*mp_rotationMatrix)*=(1+0.3*fabs(Global::msp_global->mp_normalPro->Next()));
		}else{
            mp_rotationMatrix->randomize(Program_Problem);
            mp_rotationMatrix->generateRotationMatrix(m_conditionNumber,Program_Problem);
		}
        
		ofstream out(s.c_str());
		mp_rotationMatrix->Print(out);
		out.close();
	}else{
		mp_rotationMatrix->Read_Data(in);
	}
	in.close();

    m_rotationFlag=true;
    return true;

}
double * BenchmarkFunction::getTranslation(){
    return mp_translation;
}
Matrix * BenchmarkFunction::getRotation(){
    return mp_rotationMatrix;
}
double BenchmarkFunction::getConditionNumber(){
    return m_conditionNumber;
}
ReturnFlag BenchmarkFunction::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag){
	CodeVReal &s=dynamic_cast< CodeVReal&>(ss);

	//double *x_=new double[m_numDim];
	//copy(s.m_x.begin(),s.m_x.end(),x_);
	vector<double> x_(s.m_x);
	transform(x_.data());
	evaluate__(x_.data(),s.m_obj);
	//delete [] x_;
	//x_=0;
	if (flag){
		if(rFlag)		m_evals++;	
		if(mode==Program_Algorithm&&Global::msp_global->mp_problem&&!Global::msp_global->mp_problem->isProTag(MOP)) m_globalOpt.isFound(s,m_disAccuracy,m_accuracy);
	
		if(Global::msp_global->mp_algorithm!=nullptr&&Global::msp_global->mp_algorithm->ifTerminating()){ return Return_Terminate; }
		return Return_Normal;
	}
	return Return_Normal;
} 

void BenchmarkFunction::setRotationFlag(bool rFlag){
        m_rotationFlag=rFlag;
}
void BenchmarkFunction::setTranlationFlag(bool rFlag){
    m_translationFlag=rFlag;
}
void BenchmarkFunction::setOriginalGlobalOpt(int idx,vector<double>* opt){
	m_originalGlobalOpt.setFlagLocTrue();
	if(opt==0)		for(auto&i: m_originalGlobalOpt[idx].data().m_x) i=0.;
	else	for(int i=0;i<m_numDim; i++ )  m_originalGlobalOpt[idx].data().m_x[i] = (*opt)[i];

	BenchmarkFunction::evaluate_(m_originalGlobalOpt[idx].data(),false);
}
void BenchmarkFunction::setGlobalOpt(int idx,vector<double>* opt,double *tran){
	if(m_numObj>1) throw myException("BenchmarkFunction::setGlobalOpt only for problems with a single obj");
	vector<double> v(m_numDim,0);
	m_globalOpt.setFlagLocTrue();
	if(opt==0&&tran==0)	for(int i=0;i<m_numDim;i++)v[i]+=mp_translation[i];
	else if(opt&&tran==0) for(int i=0;i<m_numDim;i++)v[i]=(*opt)[i]+mp_translation[i];
	else if(opt==0&&tran) for(int i=0;i<m_numDim;i++)v[i]+=tran[i];
	else   for(int i=0;i<m_numDim;i++)v[i]=(*opt)[i]+tran[i];
	
	m_globalOpt.setLocation(idx,v);
	BenchmarkFunction::evaluate_(m_globalOpt[idx].data(),false);
}


