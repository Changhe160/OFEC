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
// Last modified: 12 Dec. 2014
#ifndef GLOBAL_H
#define GLOBAL_H

#include "../Utility/objectFactory.h"
#include "../Problem/problem.h"

#include "../Algorithm/Algorithm.h"
#include "../Random/newran.h"
#include "../Utility/definition.h"

extern boost::mutex g_mutexStream;
extern boost::mutex g_mutex;

class Global{
public:
	unique_ptr<Cauchy> mp_cauchyPro;							// cauchy random number for algorithm
    unique_ptr<Normal> mp_normalPro;							// gaussian random number for algorithm
    unique_ptr<Levy> mp_levyPro;
    unique_ptr<Cauchy> mp_cauchyAlg;							// cauchy random number for algorithm
    unique_ptr<Normal> mp_normalAlg;							// gaussian random number for algorithm
    unique_ptr<Levy> mp_levyAlg;
	unique_ptr<Uniform> mp_uniformPro;							// random number of uniform distribution for algorithm
	unique_ptr<Uniform> mp_uniformAlg;							// random number of uniform distribution for algorithm
	unique_ptr<Problem> mp_problem;
	unique_ptr<Algorithm> mp_algorithm;						
	int m_totalNumIndis;
	int m_runId;
	
public:
	static	STRING2ID msm_pro,msm_alg;				// assign an id to each algorithm, and each problem
	static map<ALG2PRO,unsigned> msm_alg4pro;		// make pair for algorithms to problems which they solve
	#ifdef OFEC_CONSOLE
	static boost::thread_specific_ptr<Global> msp_global;
	#endif
	#ifdef OFEC_DEMON
	static unique_ptr<Global> msp_global;
	#endif
	static classFactory ms_classFactory;
	static map<string,Param> msm_param;
	static int ms_curProId,ms_curAlgId;						// Id of the problem being solved, Id of Alg runing	
	static ParamMap g_arg;
private:
	template<typename T1, typename T2>
	static bool registerItem(map<T1,T2 >& m, const T1 & name){
		if(m.size()==0 || m.size()>0&& m.end()==m.find(name) ){
			m[name]=m.size();
			return true;
		} 
		return false;
	}
public:
	Global(int runId=0,double seedPro=-1, double seedAlg=-1);
	~Global();
	void setSeedAlg(double seed);
	int getRandInt(const int min, const int max,ProgramMode mode=Program_Algorithm);
	double getRandFloat(const double min, const double max,ProgramMode mode=Program_Algorithm);

	template<typename T>
	void initializeRandomArray(T &a,const int dim,ProgramMode mode=Program_Algorithm){  // generate a set of radom numbers from 0-(dim-1) without repeat
		vector<int> temp(dim);
		for(int i=0;i<dim;i++)	temp[i]=i;
		int d=dim;
		for(int i=0;i<dim;i++){
			int t;
			if(mode==Program_Algorithm) t=	(int)(d*mp_uniformAlg->Next());
			else t=	(int)(d*mp_uniformPro->Next());
			a[i]=temp[t];
			for(int k=t;k<d-1;k++)
				temp[k]=temp[k+1];
			d--;
		}
	}
	static bool registerAlgorPro(const string&, ProgramMode flag);
	static bool registerAlg4Pro(ALG2PRO  alg2pro);
	static void registerParamter();

};

#define GET_NUM_DIM (Global::msp_global->mp_problem->getNumDim())
#define GET_NUM_OBJ (Global::msp_global->mp_problem->getNumObj())
#define GET_EVALS Global::msp_global->mp_problem->getEvaluations()

template<typename T> 
Algorithm * createAlgorithm(ParamMap &v){
	return new T(v);
}

template<typename T> 
Problem * createProblem(ParamMap &v){
	return new T(v);
}


template<typename T, ProgramMode flag>
struct RegisterClassOFEC{
	static void registerClassAlgorPro(const string &s){
		Global::ms_classFactory.m_theMapProblem.insert(make_pair(s,&createProblem<T>));
	}
};

template<typename T>
struct RegisterClassOFEC<T,Program_Algorithm>{
	static void registerClassAlgorPro(const string &s){	
		Global::ms_classFactory.m_theMapAlgorithm.insert(make_pair(s,&createAlgorithm<T>));
	}
};
template <typename T> 
int gSign(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
vector<TypeVar> gArray2Variant(const T &arr, int n){

	return vector<TypeVar>(arr,arr+n);
}

template<class T>
bool gCompare(const T d1,const T d2, bool min=true){
	if(min){
		if(d1<d2) return true;
		else return false;
	}else{
		if(d1<d2) return false;
		else return true;
	}
}

template<class T>
void gQuickSort(const T &data,int size,vector<int>& index, bool min=true,int low=0, int up=-1,int num=-1,bool start=true){
	//sort data from small to large, and put the order in index
	//size: the size of data  
	//low, up : the range of data to be sorted
	//num : the max/mim number of data within low and up 
	static boost::thread_specific_ptr<int> lb;
	static boost::thread_specific_ptr<vector<bool>> flag;
	if(start)
	{
		if(up==-1) up=size-1;
		if(num==-1) num=size;
		flag.reset(new vector<bool>(num,false));
		lb.reset(new int(low));	
		if(index.size()==0||index.size()!=size)		index.resize(size);
		for(auto i=index.begin();i!=index.end();++i) *i=i-index.begin();
	}

	
    if (low>=up) return;
	int i=0;
	for(;i<num;i++){
		if((*flag.get())[i]==false)	break;
	}
	if(i==num) return;
    int left = low+1;
    int right = up;
    int pivot=low;

    while(left<right){
		while(gCompare(data[index[left]],data[index[pivot]],min)&& left<right)         left++;
		while(!gCompare(data[index[right]],data[index[pivot]],min)&&left<right)          right--;

		 if(left<right){
            int t=index[left];
			index[left]=index[right];
			index[right]=t;
		 }
    }
	 while(!gCompare(data[index[left]],data[index[pivot]],min)&&left>pivot)  left--;
    if(gCompare(data[index[left]],data[index[pivot]],min)){
        int t=index[left];
        index[left]=index[pivot];
        index[pivot]=t;
		if(left-*lb<num)
			(*flag.get())[left-*lb]=true;
    }
	else
	{
		if(pivot-*lb<num)
			(*flag.get())[pivot-*lb]=true;
	}
	i=0;
	for(;i<num;i++){
		if((*flag.get())[i]==false)		break;
	}
	if(i==num) return;
	
    pivot=left;
	gQuickSort(data, pivot-low, index,min,low, pivot-1,num,false);
	gQuickSort(data, up-pivot,index,min, pivot+1, up,num,false);
 }
template <class T>
void gAmendSortedOrder(const T &data,int *index, int * amendedIndex,const int size){
	/*amend index in cases where same item values in data. Note: data must be
	  sorted with results stored in index, e.g.,
		data[]=[5,2,2,1]; index[]=[3,2,1,0];
		after amendation  index[]=[2,1,1,0];
	*/
	for(int r=1,idx=1;r<=size;r++,idx++) {
		int temp=r;
		int count=1;
		while(temp<size&&data[index[temp-1]]==data[index[temp]] ){count++; temp++;}
		for(int k=0;k<count;k++) amendedIndex[index[r+k-1]]=idx;
		r+=count-1;
	}
}

inline double gChaoticValue(const double x, const double min, const double max, const double rChaoticConstant=3.54){
														// return a value calculated by logistic map
	if(min>max) return -1;
	double chaotic_value;
	chaotic_value=(x-min)/(max-min);
	chaotic_value=rChaoticConstant*chaotic_value*(1-chaotic_value);
	return min+chaotic_value*(max-min);

}

string gGetProblemName(const int id);
string gGetAlgorithmName(const int id);
char *gStrtok_r(char *s, const char *delim, char **save_ptr);
bool gIsDynamicProlem();
#endif