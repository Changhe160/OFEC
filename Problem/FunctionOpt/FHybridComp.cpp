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
#include "FHybridComp.h"
#include "FSphere.h"
#include "FAckley.h"
#include "FGriewank.h"
#include "FRastrigin.h"
#include "FWeierstrass.h"
#include "FScaffer_F6.h"
#include "FElliptic.h"
#include "FNoncont_Rastrigin.h"


HybridComp::HybridComp(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){

    allocateMemory(m_numDim,m_numFuncs);
    if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05"))  {
          setSearchRange(2, 5);
		  for(int i=0;i<m_numDim; ++i)   m_searchRange[i].m_flag = false;
     }
     else  setSearchRange(-5,5);
    m_heightNormalizeSeverity=2000.;
    initialize();
}
HybridComp::HybridComp(const int rId, const int rDimNumber,  string& rName):Problem(rId,rDimNumber,rName,1),\
	BenchmarkFunction(rId,rDimNumber,rName,1){

    allocateMemory(m_numDim,m_numFuncs);
    if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05"))  {
          setSearchRange(2, 5);
		  for(int i=0;i<m_numDim; ++i)   m_searchRange[i].m_flag = false;
     }
     else  setSearchRange(-5,5);
    m_heightNormalizeSeverity=2000.;
    initialize();
}
HybridComp::~HybridComp(){
    //dtor
    freeMemory();
}
void HybridComp::allocateMemory(const int rNumDim, const int rNumFuc){
    mpp_f=new BenchmarkFunction*[rNumFuc];
    for(int i=0;i<rNumFuc;i++) mpp_f[i]=0;

    mp_convergeSeverity=new double [rNumFuc];
    mp_stretchSeverity=new double [rNumFuc];
    mp_weight=new double [rNumFuc];
    mp_height=new double[rNumFuc];
    mp_fmax=new double[rNumFuc];

}
void HybridComp::freeMemory(){
    if(m_numFuncs>0){
         for(int i=0;i<m_numFuncs;i++){
            if(mpp_f[i]) delete mpp_f[i];
         }
        delete [] mpp_f;
        delete []mp_convergeSeverity;
        delete [] mp_stretchSeverity;
        delete []mp_weight;
        delete [] mp_height;
        delete [] mp_fmax;

        mpp_f=0;
        mp_convergeSeverity=0;
        mp_stretchSeverity=0;
        mp_weight=0;
        mp_height=0;
    }


}
void HybridComp::setFunction(unsigned * rId, string rFucs[]){
    BasicFunc f;

	if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F21_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_HighConNumMatrix_F22_CEC05")\
		||IS_PROBLEM_NAME(m_id,"FUN_Noncont_RH_Com_F23_CEC05")){
		f["FUN_Expanded_Scaffer_F6_CEC05"]=&createFunction<FScaffer_F6>;
        f["FUN_Rastrigin"]=&createFunction<FRastrigin>;
        f["FUN_Griewank_Rosenbrock_F13_CEC05"]=&createFunction<FGriewank>;
        f["FUN_Weierstrass"]=&createFunction<FWeierstrass>;
        f["FUN_Griewank"]=&createFunction<FGriewank>;
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F24_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05")){
		f["FUN_Weierstrass"]=&createFunction<FWeierstrass>;
        f["FUN_Expanded_Scaffer_F6_CEC05"]=&createFunction<FScaffer_F6>;
        f["FUN_Griewank_Rosenbrock_F13_CEC05"]=&createFunction<FGriewank>;
        f["FUN_Ackley"]=&createFunction<FAckley>;
        f["FUN_Rastrigin"]=&createFunction<FRastrigin>;
        f["FUN_Griewank"]=&createFunction<FGriewank>;
        f["FUN_Noncont_Expanded_Scaffer_F6_CEC05"]=&createFunction<FScaffer_F6>;
        f["FUN_Noncont_Rastrigin"]=&createFunction<FNoncont_Rastrigin>;
        f["FUN_Elliptic"]=&createFunction<FElliptic>;
        f["FUN_Sphere_Noisy_CEC05"]=&createFunction<FSphere>;
	}else {
		f["FUN_Sphere"]=&createFunction<FSphere>;
        f["FUN_Rastrigin"]=&createFunction<FRastrigin>;
        f["FUN_Weierstrass"]=&createFunction<FWeierstrass>;
        f["FUN_Griewank"]=&createFunction<FGriewank>;
        f["FUN_Ackley"]=&createFunction<FAckley>;
            // f["Scaffer"]=&createFunction<FScaffer>;
	}

    for(int i=0;i<m_numFuncs;i++){
		mpp_f[i]=dynamic_cast<BenchmarkFunction*>(f[rFucs[i]](rId[i],m_numDim,rFucs[i]));
        mpp_f[i]->setBias(0);
    }

	if(IS_PROBLEM_NAME(m_id,"FUN_R_Com")||IS_PROBLEM_NAME(m_id,"FUN_H_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_H_Com_Noisy_CEC05")||\
		IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F21_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_Noncont_RH_Com_F23_CEC05")){
		for(int i=0;i<m_numFuncs;i++) mpp_f[i]->setConditionNumber(2.);
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_HighConNumMatrix_F22_CEC05")){
		mpp_f[0]->setConditionNumber(10); mpp_f[1]->setConditionNumber(20);
        mpp_f[2]->setConditionNumber(50);mpp_f[3]->setConditionNumber(100);
        mpp_f[4]->setConditionNumber(200);mpp_f[5]->setConditionNumber(1000);
        mpp_f[6]->setConditionNumber(2000);mpp_f[7]->setConditionNumber(3000);
        mpp_f[8]->setConditionNumber(4000);mpp_f[9]->setConditionNumber(5000);
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_NarrowBasin_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_Bound_CEC05")){
		mpp_f[0]->setConditionNumber(2); mpp_f[1]->setConditionNumber(3);
        mpp_f[2]->setConditionNumber(2);mpp_f[3]->setConditionNumber(3);
        mpp_f[4]->setConditionNumber(2);mpp_f[5]->setConditionNumber(3);
        mpp_f[6]->setConditionNumber(20);mpp_f[7]->setConditionNumber(30);
        mpp_f[8]->setConditionNumber(200);mpp_f[9]->setConditionNumber(300);
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F24_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05")){
		mpp_f[0]->setConditionNumber(100); mpp_f[1]->setConditionNumber(50);
		mpp_f[2]->setConditionNumber(30);mpp_f[3]->setConditionNumber(10);
		mpp_f[4]->setConditionNumber(5);mpp_f[5]->setConditionNumber(5);
		mpp_f[6]->setConditionNumber(4);mpp_f[7]->setConditionNumber(3);
		mpp_f[8]->setConditionNumber(2);mpp_f[9]->setConditionNumber(2);
	}else{
		for(int i=0;i<m_numFuncs;i++) mpp_f[i]->setConditionNumber(2.);
	}

}

int HybridComp::getNumFuncs(){
    return m_numFuncs;
}

bool HybridComp::loadTranslation(){
	string s;
	stringstream ss;
	ss<<m_numDim<<"Dim.txt";
	s=ss.str();
	s.insert(0,m_name+"_Opt_");

	s.insert(0,"Problem/FunctionOpt/Data/");//probDataPath
	s.insert(0,Global::g_arg[param_workingDir]);
	ifstream in;
	in.open(s.c_str());

	if(in.fail()){
		if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_NarrowBasin_CEC05")){
			for(int i=0;i<m_numFuncs-1;i++)
			for(int j=0;j<m_numDim;j++){
				mpp_f[i]->getTranslation()[j]=m_searchRange[j].m_lower +(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*(1-Global::msp_global->mp_uniformPro->Next());
			}

			for(int j=0;j<m_numDim;j++){
				mpp_f[m_numFuncs-1]->getTranslation()[j]=0;
			}
		}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_Bound_CEC05")){
			for(int i=0;i<m_numFuncs-1;i++)
			for(int j=0;j<m_numDim;j++){
				if(i==0&&(j+1)%2==0) mpp_f[i]->getTranslation()[j]=m_searchRange[j].m_upper;
				else mpp_f[i]->getTranslation()[j]=m_searchRange[j].m_lower +(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*(1-Global::msp_global->mp_uniformPro->Next());
			}

			for(int j=0;j<m_numDim;j++){
				mpp_f[m_numFuncs-1]->getTranslation()[j]=0;
			}
		}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05")){
			for(int j=0;j<m_numDim;j++){
				mpp_f[0]->getTranslation()[j]= -5.0 +(7.0)*(Global::msp_global->mp_uniformPro->Next());
               }

			for(int i=1;i<m_numFuncs-1;i++)
			for(int j=0;j<m_numDim;j++){
				mpp_f[i]->getTranslation()[j]=m_searchRange[j].m_lower +(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*(1-Global::msp_global->mp_uniformPro->Next());
			}

			for(int j=0;j<m_numDim;j++){
				mpp_f[m_numFuncs-1]->getTranslation()[j]=0;
			}
		}else{
			for(int i=0;i<m_numFuncs;i++)
			for(int j=0;j<m_numDim;j++){
				mpp_f[i]->getTranslation()[j]=m_searchRange[j].m_lower +(m_searchRange[j].m_upper-m_searchRange[j].m_lower)*(1-Global::msp_global->mp_uniformPro->Next());
            }
		}

		
		ofstream out(s.c_str());
		for(int i=0;i<m_numFuncs;i++){
			for(int j=0;j<m_numDim;j++)
				out<<mpp_f[i]->getTranslation()[j]<<" ";
				out<<endl;
		}
		out.close();

	}else{
        for(int i=0;i<m_numFuncs;i++){
			for(int j=0;j<m_numDim;j++) {
				in>>mpp_f[i]->getTranslation()[j];
			}
		}
	}
	in.close();

    for(int i=0;i<m_numFuncs;i++)	mpp_f[i]->setTranlationFlag(true);

    return true;
}

bool HybridComp::loadRotation(){
    string s;
	char a[100];
	sprintf(a,"%d",m_numDim);
	strcat(a,"Dim.txt");
	s=a;
	s.insert(0,m_name+"_RotM_");

	s.insert(0,"Problem/FunctionOpt/Data/");//probDataPath
	s.insert(0,Global::g_arg[param_workingDir]);
	ifstream in;
	in.open(s.c_str());
	if(in.fail()){

		for(int i=0;i<m_numFuncs;i++){

			if(IS_PROBLEM_NAME(m_id,"FUN_Com")||IS_PROBLEM_NAME(m_id,"FUN_Com_CEC05")){
				 mpp_f[i]->getRotation()->identity();
			}else{
				mpp_f[i]->getRotation()->randomize(Program_Problem);
				mpp_f[i]->getRotation()->generateRotationMatrix(mpp_f[i]->getConditionNumber(),Program_Problem);
			}
		}
		ofstream out(s.c_str());
		for(int i=0;i<m_numFuncs;i++) mpp_f[i]->getRotation()->Print(out);
		out.close();
	}else{
	    for(int i=0;i<m_numFuncs;i++)	mpp_f[i]->getRotation()->Read_Data(in);
	}

	for(int i=0;i<m_numFuncs;i++)	mpp_f[i]->setRotationFlag(true);
	in.close();
    return true;
}
void HybridComp::initialize(){

	 m_originalGlobalOpt.flagGloObj()=false;
	 m_originalGlobalOpt.flagLoc()=false;

	if(IS_PROBLEM_NAME(m_id,"FUN_Com")||IS_PROBLEM_NAME(m_id,"FUN_R_Com")){
		setUpFCom();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_H_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_H_Com_Noisy_CEC05")){
		SetUpFCom_CEC05();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_NarrowBasin_CEC05")||\
		IS_PROBLEM_NAME(m_id,"FUN_RH_Com_Bound_CEC05")){
		SetUpFRH_Com_CEC05();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F21_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_HighConNumMatrix_F22_CEC05")||\
		IS_PROBLEM_NAME(m_id,"FUN_Noncont_RH_Com_F23_CEC05")){
		SetUpFRH001_Com_CEC05();
	}else if(IS_PROBLEM_NAME(m_id,"FUN_RH_Com_F24_CEC05")||IS_PROBLEM_NAME(m_id,"FUN_RH_Com_noBounds_F25_CEC05")){
		 SetUpFRH002_Com_CEC05();
	}else {
		 throw myException("Error: please check the problem ID@HybridComp::initialize");
	}

    loadTranslation();
    loadRotation();

	CodeVReal x(m_numDim,m_numObj);
    for(int i=0;i<m_numFuncs;i++){
        for(int j=0;j<m_numDim;j++){ // calculate the estimate max value of funciton i
            x[j]=m_searchRange[j].m_upper;
            x[j]/=mp_stretchSeverity[i];
        }		
		mpp_f[i]->BenchmarkFunction::evaluate_(x,false);
		mp_fmax[i]=x.m_obj[0];
     }
    
	vector<double> v(m_numDim,0);
	m_globalOpt.setNumOpts(1);
    setGlobalOpt(0,&v,mpp_f[0]->getTranslation());
    setAccuracy(1.0e-3);
	setDisAccuracy(0.1);
}

void HybridComp::SetUpFRH002_Com_CEC05()
{
     string fname[m_numFuncs];
    unsigned fid[m_numFuncs];
	for(int i=0;i<m_numFuncs;i++){
		mp_height[i]=100*i;
	}
	fname[0] = "FUN_Weierstrass";
	fname[1] = "FUN_Expanded_Scaffer_F6_CEC05";
	fname[2] = "FUN_Griewank_Rosenbrock_F13_CEC05";
	fname[3] = "FUN_Ackley";
	fname[4] = "FUN_Rastrigin";
	fname[5] = "FUN_Griewank";
	fname[6] = "FUN_Noncont_Expanded_Scaffer_F6_CEC05";
	fname[7] = "FUN_Noncont_Rastrigin";
	fname[8] = "FUN_Elliptic";
	fname[9] = "FUN_Sphere_Noisy_CEC05";

	mp_stretchSeverity[0]=10;
     mp_convergeSeverity[0]=2.;

	 for(int i=0;i<m_numFuncs;i++)   fid[i]=Global::msm_pro[fname[i]];
    setFunction(fid,fname);

    mpp_f[0]->setSearchRange(-0.5,0.5);
    mpp_f[1]->setSearchRange(-100, 100);
	mpp_f[2]->setSearchRange(-5,5);
	mpp_f[3]->setSearchRange(-32,32);
	mpp_f[4]->setSearchRange(-5,5);
	mpp_f[5]->setSearchRange(-5,5);
    mpp_f[6]->setSearchRange(-100, 100);
    mpp_f[7]->setSearchRange(-5,5);
	mpp_f[8]->setSearchRange(-100,100);
	mpp_f[9]->setSearchRange(-100,100);

	mp_stretchSeverity[1]=5.0/20;
	mp_stretchSeverity[2]=1.0;
	mp_stretchSeverity[3]=5.0/32;
	mp_stretchSeverity[4]=1.0;
	mp_stretchSeverity[5]=5.0/100;
	mp_stretchSeverity[6]=5.0/50;
	mp_stretchSeverity[7]=1.0;
	mp_stretchSeverity[8]=5.0/100;
	mp_stretchSeverity[9]=5.0/100;

	mp_convergeSeverity[1]=2.;
	mp_convergeSeverity[2]=2.;	mp_convergeSeverity[3]=2.;
	mp_convergeSeverity[4]=2.;  mp_convergeSeverity[5]=2.;
	mp_convergeSeverity[6]=2;	mp_convergeSeverity[7]=2;
	mp_convergeSeverity[8]=2.;  mp_convergeSeverity[9]=2.;

    for(int i=0;i<m_numFuncs;i++){
		mpp_f[i]->setScale(mp_stretchSeverity[i]);
	}
	setBias(260);
}

void HybridComp::setUpFCom(){

	for(int i=0;i<m_numFuncs;i++){
		mp_height[i]=100*i;
		mp_convergeSeverity[i]=1.;
	}

	string fname[m_numFuncs];
    unsigned fid[m_numFuncs];


	fname[0]="FUN_Sphere";		fname[1]="FUN_Sphere";
	fname[2]="FUN_Rastrigin";	fname[3]="FUN_Rastrigin";
	fname[4]="FUN_Weierstrass";	fname[5]="FUN_Weierstrass";
	fname[6]="FUN_Griewank";	fname[7]="FUN_Griewank";
	fname[8]="FUN_Ackley";		fname[9]="FUN_Ackley";

    for(int i=0;i<m_numFuncs;i++)  {
		fid[i]=Global::msm_pro[fname[i]];
    }

    setFunction(fid,fname);

	mpp_f[0]->setSearchRange(-100,100); mpp_f[1]->setSearchRange(-100,100);
	mpp_f[2]->setSearchRange(-5,5);     mpp_f[3]->setSearchRange(-5,5);
	mpp_f[4]->setSearchRange(-0.5,0.5); mpp_f[5]->setSearchRange(-0.5,0.5);
	mpp_f[6]->setSearchRange(-100,100); mpp_f[7]->setSearchRange(-100,100);
	mpp_f[8]->setSearchRange(-32,32);   mpp_f[9]->setSearchRange(-32,32);


	for(int i=0;i<m_numFuncs;i++){
		double l,u;
		mpp_f[i]->getSearchRange(l,u,0);
		mp_stretchSeverity[i]=mp_convergeSeverity[i]*(m_searchRange[0].m_upper-m_searchRange[0].m_lower)/(u-l);
		mpp_f[i]->setScale(mp_stretchSeverity[i]);
	}
   // setBias(0);
}

void HybridComp::SetUpFRH001_Com_CEC05()
{
     string fname[m_numFuncs];
	 unsigned fid[m_numFuncs];
	for(int i=0;i<m_numFuncs;i++){
		mp_height[i]=100*i;
	}
	fname[0]="FUN_Expanded_Scaffer_F6_CEC05";			fname[1]="FUN_Expanded_Scaffer_F6_CEC05";
	fname[2]="FUN_Rastrigin";							fname[3]="FUN_Rastrigin";
	fname[4]="FUN_Griewank_Rosenbrock_F13_CEC05";		fname[5]="FUN_Griewank_Rosenbrock_F13_CEC05";
	fname[6]="FUN_Weierstrass";							fname[7]="FUN_Weierstrass";
	fname[8]="FUN_Griewank";							fname[9]="FUN_Griewank";

	mp_stretchSeverity[0]=5.*5./100;
     mp_convergeSeverity[0]=1.;

    for(int i=0;i<m_numFuncs;i++)   fid[i]=Global::msm_pro[fname[i]];
    setFunction(fid,fname);
    mpp_f[0]->setSearchRange(-100,100);   mpp_f[1]->setSearchRange(-100,100);
	mpp_f[2]->setSearchRange(-5,5);     mpp_f[3]->setSearchRange(-5,5);
	mpp_f[4]->setSearchRange(-5,5); mpp_f[5]->setSearchRange(-5,5);
    mpp_f[6]->setSearchRange(-0.5,0.5); mpp_f[7]->setSearchRange(-0.5,0.5);
	mpp_f[8]->setSearchRange(-200,200); mpp_f[9]->setSearchRange(-200,200);

	mp_stretchSeverity[1]=5./100;
	mp_stretchSeverity[2]=5.;		mp_stretchSeverity[3]=1.;
	mp_stretchSeverity[4]=5;  mp_stretchSeverity[5]=1;
	mp_stretchSeverity[6]=50.;	mp_stretchSeverity[7]=10.;
	mp_stretchSeverity[8]=5.*5/200;  mp_stretchSeverity[9]=5./200;

	mp_convergeSeverity[1]=1.;
	mp_convergeSeverity[2]=1.;	mp_convergeSeverity[3]=1.;
	mp_convergeSeverity[4]=1.;  mp_convergeSeverity[5]=2.;
	mp_convergeSeverity[6]=2;	mp_convergeSeverity[7]=2;
	mp_convergeSeverity[8]=2.;  mp_convergeSeverity[9]=2.;

    for(int i=0;i<m_numFuncs;i++){
		mpp_f[i]->setScale(mp_stretchSeverity[i]);
	}
	setBias(360);
}

void HybridComp::SetUpFCom_CEC05(){
	string fname[m_numFuncs];
    unsigned fid[m_numFuncs];
	for(int i=0;i<m_numFuncs;i++){
		mp_height[i]=100*i;
		mp_convergeSeverity[i]=1.;
	}
	fname[0]="FUN_Rastrigin";	fname[1]="FUN_Rastrigin";
	fname[2]="FUN_Weierstrass";	fname[3]="FUN_Weierstrass";
	fname[4]="FUN_Griewank";	fname[5]="FUN_Griewank";
	fname[6]="FUN_Ackley";		fname[7]="FUN_Ackley";
	fname[8]="FUN_Sphere";		fname[9]="FUN_Sphere";

    for(int i=0;i<m_numFuncs;i++)   {
         fid[i]=Global::msm_pro[fname[i]];
    }
    setFunction(fid,fname);


	mpp_f[0]->setSearchRange(-5,5);     mpp_f[1]->setSearchRange(-5,5);
	mpp_f[2]->setSearchRange(-0.5,0.5); mpp_f[3]->setSearchRange(-0.5,0.5);
	mpp_f[4]->setSearchRange(-60,60); mpp_f[5]->setSearchRange(-60,60);
	mpp_f[6]->setSearchRange(-32,32);   mpp_f[7]->setSearchRange(-32,32);
	mpp_f[8]->setSearchRange(-100,100); mpp_f[9]->setSearchRange(-100,100);


	mp_stretchSeverity[0]=1.;		mp_stretchSeverity[1]=1.;
	mp_stretchSeverity[2]=10.;		mp_stretchSeverity[3]=10.;
	mp_stretchSeverity[4]=5./60;  mp_stretchSeverity[5]=5./60;
	mp_stretchSeverity[6]=5./32;	mp_stretchSeverity[7]=5./32;
	mp_stretchSeverity[8]=5./100;  mp_stretchSeverity[9]=5./100;

    for(int i=0;i<m_numFuncs;i++){
		mpp_f[i]->setScale(mp_stretchSeverity[i]);
	}

    setBias(120.);
}
void HybridComp::SetUpFRH_Com_CEC05(){
	string fname[m_numFuncs];
	unsigned fid[m_numFuncs];
	for(int i=0;i<m_numFuncs;i++){
		mp_height[i]=100*i;

	}
	fname[0]="FUN_Ackley";			fname[1]="FUN_Ackley";
	fname[2]="FUN_Rastrigin";		fname[3]="FUN_Rastrigin";
	fname[4]="FUN_Sphere";			fname[5]="FUN_Sphere";
	fname[6]="FUN_Weierstrass";		fname[7]="FUN_Weierstrass";
	fname[8]="FUN_Griewank";		fname[9]="FUN_Griewank";

	if(IS_PROBLEM_NAME( m_id,"FUN_RH_Com_NarrowBasin_CEC05"))   {
	    mp_stretchSeverity[0]=0.1*5./32;
		mp_convergeSeverity[0]=0.1;
	}else{
        mp_stretchSeverity[0]=2.*5./32;
		mp_convergeSeverity[0]=1.;
	}
	for(int i=0;i<m_numFuncs;i++)   fid[i]=Global::msm_pro[fname[i]];
    setFunction(fid,fname);


    mpp_f[0]->setSearchRange(-32,32);   mpp_f[1]->setSearchRange(-32,32);
	mpp_f[2]->setSearchRange(-5,5);     mpp_f[3]->setSearchRange(-5,5);
	mpp_f[4]->setSearchRange(-100,100); mpp_f[5]->setSearchRange(-100,100);
    mpp_f[6]->setSearchRange(-0.5,0.5); mpp_f[7]->setSearchRange(-0.5,0.5);
	mpp_f[8]->setSearchRange(-60,60); mpp_f[9]->setSearchRange(-60,60);

	mp_stretchSeverity[1]=5./32;
	mp_stretchSeverity[2]=2.;		mp_stretchSeverity[3]=1.;
	mp_stretchSeverity[4]=2*5./100;  mp_stretchSeverity[5]=5./100;
	mp_stretchSeverity[6]=20.;	mp_stretchSeverity[7]=10.;
	mp_stretchSeverity[8]=2.*5/60;  mp_stretchSeverity[9]=5./60;

	mp_convergeSeverity[1]=2.;
	mp_convergeSeverity[2]=1.5;	mp_convergeSeverity[3]=1.5;
	mp_convergeSeverity[4]=1.;  mp_convergeSeverity[5]=1.;
	mp_convergeSeverity[6]=1.5;	mp_convergeSeverity[7]=1.5;
	mp_convergeSeverity[8]=2.;  mp_convergeSeverity[9]=2.;

    for(int i=0;i<m_numFuncs;i++){
		mpp_f[i]->setScale(mp_stretchSeverity[i]);
	}
	setBias(10.);
}

 void HybridComp::evaluate__(double  const *x,vector<double>& obj_){
     double * x_ = new double[m_numDim];
      if(IS_PROBLEM_NAME(m_id,"FUN_Noncont_RH_Com_F23_CEC05")) {
           for(int j=0;j<m_numDim; ++j) {
               double trans = mpp_f[0]->getTranslation()[j];
               if(fabs(x[j] - trans) >= 0.5) {
                    double xTemp = 2.0 * x[j];
                    int intPart = int(xTemp);
                    double decimalPart = xTemp - intPart;
                    if(xTemp <= 0 && decimalPart >= 0.5) xTemp = intPart - 1;
                    else if(decimalPart < 0.5) xTemp = intPart;
                    else xTemp = intPart + 1;
                    x_[j] = xTemp / 2.0;
               }
               else x_[j] = x[j];
           }
      }
	  else copy(x,x+m_numDim,x_); 
    vector<double> weight(m_numFuncs,0);
	for(int i=0;i<m_numFuncs;i++){ // calculate mp_weight for each function
		for(int j=0;j<m_numDim;j++) {
			weight[i]+=(x_[j]-mpp_f[i]->getTranslation()[j])*(x_[j]-mpp_f[i]->getTranslation()[j]);
		}
		weight[i]=exp(-sqrt(weight[i]/(2*m_numDim*mp_convergeSeverity[i]*mp_convergeSeverity[i])));
	}
	vector<double> fit(m_numFuncs);
	CodeVReal s(m_numDim,m_numObj);
	for(int i=0;i<m_numFuncs;i++){ // calculate objective value for each function
		copy(x_,x_+m_numDim,s.m_x.begin());
		mpp_f[i]->BenchmarkFunction::evaluate_(s,false);
		fit[i]=s.m_obj[0];
		if(mp_fmax[i]!=0)
		fit[i]=m_heightNormalizeSeverity*fit[i]/fabs(mp_fmax[i]);
	}

	double sumw=0,wmax;
	wmax=*max_element(weight.begin(),weight.end());
	for(int i=0;i<m_numFuncs;i++) {
		if(weight[i]!=wmax) {
			weight[i]=weight[i]*(1-pow(wmax,10));
		}
	}
     int sameWmax_N = 0;
     for(int i=0; i<m_numFuncs; ++i) {
          if(weight[i] == wmax) ++sameWmax_N;
     }
     int i=m_numFuncs-1;
     while(sameWmax_N > 1 && i >= 0) {
          if(wmax == weight[i]) {
               weight[i] = 0;
               --sameWmax_N;
          }
          --i;
     }

	for(int i=0;i<m_numFuncs;i++)
		sumw+=weight[i];
	for(int i=0;i<m_numFuncs;i++)
		weight[i]/=sumw;

	double obj=0;
	for(int i=0;i<m_numFuncs;i++) {
		obj+=weight[i]*(fit[i]+mp_height[i]);
	}

	delete[] x_;
	x_=0;

	if(IS_PROBLEM_NAME(m_id,"FUN_H_Com_Noisy_CEC05")) obj_[0]= obj*(1+0.2*fabs(Global::msp_global->mp_normalPro->Next()))+m_bias;
	else	obj_[0]=obj+m_bias;

 }



