/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Oct 2014
// Last modified:

#ifndef  TRAVELLING_SALESMAN_H
#define  TRAVELLING_SALESMAN_H

#include "../../problem.h"
#include "../../optimum.h"
#include "../../../Global/boundary.h"
#define CAST_TSP dynamic_cast<TravellingSalesman *>(Global::msp_global->mp_problem.get())
class TravellingSalesman :public Problem
{
protected:
	vector<vector< vector<double> > >  mvvv_cost;       //the cost between city and city
	string m_fileName;
	vector<vector<bool>> mvv_solut;
	BoundaryTSP m_searchRange;
	Optima<CodeVInt> m_globalOpt; 
public:
	TravellingSalesman(ParamMap& v);
	~TravellingSalesman();
	TravellingSalesman(const int rId, const int rDimNumber, string rName, string fileName, int numObj=1);
	ReturnFlag evaluate_(VirtualEncoding &s, bool rFlag, ProgramMode mode = Program_Problem, bool flag=true);
	bool isValid(const VirtualEncoding &s);
	void initializeSolution(VirtualEncoding &result,const int idx=0,const int maxId=0);
	void initializeSolution(const VirtualEncoding &base,VirtualEncoding &result,double range){}
	void initializePartSolution(VirtualEncoding &result,int begin,int end){}
	bool isSame(const VirtualEncoding &s1, const VirtualEncoding &s2);

	virtual void readProblem();    //read source data from file
	void findNearbyCity(vector<vector<int> > &nearby,int n=0);     //find some percent of nearby city
	virtual void createCandidateSets(vector<vector<int> > &candidateSets);  //create candidate sets
	void prim(vector<vector<int> > &mstEdge,int n=0); //find MST edges
	void calculateEdgeWeight(char *edgeType,vector<vector<double> >& coordinate);
	const vector< vector<double> > & getCost(int i=0) const { return mvvv_cost[i]; }
	string getFileName() const { return m_fileName; }
	double getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode);
	pair<int,int> getNextCity(const VirtualEncoding &s,int n);  //return the next node of node n in solution s

	bool getObjGlobalOpt(vector<double> &opt);
	bool getObjGlobalOpt(vector<vector<double>> &opt);
	const Optima<CodeVInt> & getGOpt()const; 
	Optima<CodeVInt> & getGOpt();
	bool isGlobalOptKnown();
	TravellingSalesman *getTypePtr();
	TravellingSalesman &getTypeRef();
	void validate(VirtualEncoding &s,SolutionValidation *mode=0){}
};

int selectCityRandom(vector< vector<int> > &,vector<int> &visited,int num,int row);
int selectCityRandom(vector<int> &visited,int dim);

#endif  //TRAVELLING_SALESMAN_H