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
// Created: 9 Aug 2014
// Last modified:
#ifndef PROBLEM_H
#define PROBLEM_H

#include "../Global/encoding.h"

extern int &g_curProID;
class Problem{
    public:
        string m_name;
        stringstream m_proPar;
     protected:
        const int m_id;
		int m_numDim;
		vector<Compare> m_OptMode;
        int m_evals;
		int m_numObj;
		double m_accuracy;
		bool m_sameType;
		int m_tevals,m_cevals;	// the total evals and the number of countable evlas 
		
		PopInitMethod m_popInitialMode;
		SolutionValidation m_validationMode;
		set<ProTag> m_tag;
		vector<pair<double, double>>m_objRange;
		vector<vector<double>*> m_os;  //the set of objectives
	protected:
		virtual void  freeMemory(){};
		virtual void parameterSetting(Problem * rDP);
		void setSameType(bool flag);
		Problem():m_id(-1){}
		void setObjNumber(const int ObjNum);	
		virtual void allocateMemory(const int numDim){}
		virtual void resizeDim(int num){}
		virtual void resizeObj(int num);
		virtual void setObjSet() = 0;
	public:
		Problem(const int rId, const int rDimNumber, string rName, const int numObj=1);
		virtual ~Problem(){};
		Problem& operator=(const Problem & rP);
		int getId() const;
		inline int getNumDim() const;
		Compare getOptType(int idx=0) const;
		inline int getEvaluations() const;
		void resetEvaluations();

		void setOptType(const Compare rT, int idx = 0);
		inline int getNumObj();
		virtual ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true) = 0;
		template<typename SOL>
		ReturnFlag evaluate(SOL &s, bool rFlag=true,ProgramMode mode=Program_Algorithm ); 	
		
		virtual void reset(){}
		virtual void reinitialize(){}

		virtual bool isValid(const VirtualEncoding &s)=0;
		virtual void validate(VirtualEncoding &s,SolutionValidation *mode=0)=0;
		virtual void initializeSolution(VirtualEncoding &result,const int idx=0,const int maxId=0)=0;
		virtual void initializeSolution(const VirtualEncoding &ref,VirtualEncoding &result,double range)=0;
		virtual void initializePartSolution(VirtualEncoding &result,int begin,int end)=0;
		void setPopInitialMode(PopInitMethod m);
		void setValidateMode(SolutionValidation m);
		virtual double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode)=0;
		virtual bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2)=0;
		virtual bool isGlobalOptKnown()=0;
		inline double getAccuracy();
		void setAccuracy(double rAcc);
		bool isSameType();

		void setProTag(const set<ProTag> &);
		void addProTag(ProTag tag);
		bool isProTag(ProTag tag);
		inline int &cevals();
		inline int getTevals();

		virtual bool getObjGlobalOpt(vector<double> &opt)=0;
		virtual bool getObjGlobalOpt(vector<vector<double>> &opt)=0;
		virtual CompareResultFlag compare(const VirtualEncoding &s1, const VirtualEncoding &s2)const;
		virtual void copyChanges(const Problem * pro, const vector<int> *cd = nullptr, const vector<int> *co = nullptr);
		virtual const vector<pair<double, double>>& getObjRange() { return m_objRange; }
		const vector<vector<double>*> & getObjSet() { return m_os; }

		virtual bool isGlobalOptFound() = 0;
};

template<typename T> Problem * createFunction( int rId,  int rDimNumber, string& rName){
	return new T( rId,  rDimNumber, rName);
}
typedef Problem *(*pFun)(int rId,  int rDimNumber, string& rName);
typedef map<string, pFun > BasicFunc;


template<typename SOL>
ReturnFlag  Problem::evaluate(SOL &s, bool rFlag,ProgramMode mode ){
	bool inbound=true;
	if(m_id==g_curProID/*Global::ms_curProId*/){
		inbound=isValid(s.data());
    }
	
	if (!inbound){
		#ifdef OFEC_DEMON
		if (m_evals == 0){		
			for (int i = 0; i < m_numObj;++i) s.obj()[i] = (m_OptMode[i] == MIN_OPT) ? LONG_MAX : LONG_MIN;
		}
		#endif
		#ifdef OFEC_CONSOLE
		for (int i = 0; i < m_numObj;++i) s.obj()[i] = (m_OptMode[i] == MIN_OPT) ? LONG_MAX : LONG_MIN;
		#endif
		return Return_Error;
	}

	ReturnFlag rf=evaluate_(s.data(),rFlag,mode);
	if(rFlag){
		m_cevals++;
	}
	++m_tevals;


	return rf;
	
}

int & Problem::cevals(){
	return m_cevals;
}

int  Problem::getTevals(){
	return m_tevals;
}

double Problem::getAccuracy(){ 
	return m_accuracy; 
}

int Problem::getEvaluations() const{      
	return m_evals;       
}

int Problem::getNumDim() const{   
	return m_numDim;       
}

int Problem::getNumObj(){ 
	return m_numObj;
}
#endif // PROBLEM_H
