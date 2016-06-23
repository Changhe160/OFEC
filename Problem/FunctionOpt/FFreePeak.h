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
*******************************************************************************************/
// Created: 21 June 2015
// Last modified:
#ifndef FFREEPEAK_H
#define FFREEPEAK_H
#include "BenchmarkFunction.h"
#include "../../Utility/kdtree_space.h"
class FreePeak :public ContinuousProblem
{
public:
	enum DivisionMode{DM_Random=1,DM_Even,DM_UserDefined};
	struct Box{
		vector < unique_ptr<BenchmarkFunction> > peak;	//each peak must be a single-obj unimodal problem
		long double boxSize;
		double boxRatio;
	};

	FreePeak(ParamMap &v);
	~FreePeak(){};
	FreePeak(const FreePeak&) = delete;
	const CodeVReal & getPeak(int idx) {
		return m_peak[idx];
	}
	const CodeVReal & getPeak(int bidx,int oidx) {
		return m_peak[m_peakIdx[make_pair(bidx,oidx)]];
	}
	const vector<CodeVReal> &getPeak() {
		return m_peak;
	}
	FreePeak& operator = (FreePeak& rhs);
	void copyChanges(const Problem * pro, const vector<int> *cd = nullptr, const vector<int> *co = nullptr);
	const vector < unique_ptr<BenchmarkFunction> >& getBox(int idx){ return (*m_box)[idx].peak; }
	int getNumBox(){ return m_numBox; }
	vector<double> mapToPeak(const vector<double> &p, int gidx = -1, const int pidx = 0, const int tidx = 0);
	vector<double> mapFromPeak(const vector<double> &p, const int gidx, const int pidx = 0, const int tidx=0);
	void setBoxSize(const vector<double>& rat);
	void initializeSolution(VirtualEncoding &result, const int idx = 0, const int maxId = 0);
	void setInitialBox4Solution(const vector<int>&v){ m_initialBox4Sol = v; };
	int getNumPeak(){ return m_numPeak; }
	KDTreeSpace::PartitioningKDTree<double> *getTree(int idx=0){ return m_tree[idx].get(); }
	int getNumTree(){ return m_tree.size(); }
	int getTreeRoot(const vector<double> &p)const;
	int boxIdx(const int tidx, const int gidx)const;
	void treeIdx(const int bidx, int &tidx, int &gidx)const;
	int treeIdx(const int bidx)const;

protected:	
	unique_ptr<vector<Box>> m_box;
	int m_numPeak=0;
	int m_numBox;
	vector<vector<vector<double>>> m_divisionPoint;
	vector<unique_ptr<KDTreeSpace::PartitioningKDTree<double>>> m_tree;
	int m_peaksPerBox=1;
	vector<BenchmarkFunction*> m_fp;
	map<tuple<int, int,int>, int> m_objIdx;
	map<pair<int, int>, int> m_peakIdx;
	vector<CodeVReal> m_peak;
	int m_dm = DivisionMode::DM_Random;
	vector<int>m_initialBox4Sol;
protected:
	void initialize();
	void createPeaks();
	virtual void createPeaks_() = 0;
	double disNearestPeak();
	void computeBoxSize();
	void updateGOpt();
	void updatePeak();	
	void resizeDim(int num);
	void resizeObj(int num);

	int smallestBox(const set<int>& except);
	int largestBox(const set<int>& except);

public:
	struct Peak {
		int shape;
		string basin;
		vector<int> transf;
		double height, minHeigh;
	};
	static void generateDivision(const vector<vector<pair<double, double>>> &box, vector<int> &node, int numNode, int numDim, int numTree, const string &path, bool evendiv=false,vector<double> *br=0);
	static void generateLocationSingleObj(int numDim, int numPeak, const string & path, bool random=true);  // mode=true denotes that peak locations are randomly generated in the space, otherwise at the center
	static void generateLocationMultiObj(int numDim, int numPeak, int numObj,float minRadius /*[0,1]*/,const string & path,int mode, float rhoweb=1.); //mode=1: Jump, mode=2: Web 
	static void generateConf(const vector<Peak> &conf, const vector<int> &gbox, const string & path, int obj = 1);
	static void generateConfGOP(int numPeak, const string & path, int gShape, int gTrans,const string & gPosition,int lShape,int lTrans,double lMaxH, bool lRandH, bool trap, double tMinH,double tMaxH);
	static void generateConfMMOP(int numPeak, const string & path, int gShape, int gTrans, const string & gPosition, int gNum, int lShape, int lTrans, double lMaxH, bool lRandH, bool trap, int tNum,double tMinH, double tMaxH);
	static void generateConfMOP(int numBox, const string & path, int o1stShape, int numObj, int mode, float rhoHeight, int numParetoRegion = 1); //mode=1: Jump, mode=2: Web, mode=3: Countable
	static void setup(const string & problem);
};

//p is a point in free peak
inline vector<double> FreePeak::mapToPeak(const vector<double> &p, int gidx, const int pidx, const int tidx ){
	if (gidx == -1) gidx = m_tree[tidx]->get_regionIdx(p);
	int idx = boxIdx(tidx, gidx);
	vector<double> r(m_numDim);
	double l, u;
	for (int i = 0; i < m_numDim; ++i){
		(*m_box)[idx].peak[pidx]->getSearchRange(l, u, i);
		r[i] = l + (p[i] - m_tree[tidx]->region[gidx].box[i].first) / (m_tree[tidx]->region[gidx].box[i].second - m_tree[tidx]->region[gidx].box[i].first)*(u - l);
	}
	return move(r);
}
//p is a point in a sub-problem
inline vector<double> FreePeak::mapFromPeak(const vector<double> &p, const int gidx, const int pidx, const int tidx){
	vector<double> r(m_numDim);
	int idx = boxIdx(tidx, gidx);
	double l, u;
	for (int i = 0; i < m_numDim; ++i){
		(*m_box)[idx].peak[pidx]->getSearchRange(l, u, i);
		r[i] = m_tree[tidx]->region[gidx].box[i].first + (p[i] - l) / (u - l)*(m_tree[tidx]->region[gidx].box[i].second - m_tree[tidx]->region[gidx].box[i].first);
	}
	return move(r);
}

inline int FreePeak::getTreeRoot(const vector<double> &p)const{
	if (m_tree.size() == 1) return 0;

	for (int i = 0; i < m_tree.size(); ++i){
		bool flag = true;
		for (int j = 0; j < m_numDim; ++j){
			if (m_tree[i]->get_rootBox().box[j].second == m_searchRange[j].m_upper){
				if (p[j] < m_tree[i]->get_rootBox().box[j].first){
					flag = false;
					break;
				}
			}
			else{
				if (p[j] < m_tree[i]->get_rootBox().box[j].first || p[j] >= m_tree[i]->get_rootBox().box[j].second){
					flag = false;
					break;
				}
			}			
		}
		if (flag) return i;
	}

	return -1;
}

inline int FreePeak::boxIdx(const int tidx, const int gidx)const{
	int idx = 0;
	for (int i = 0; i < tidx; ++i){
		idx += m_divisionPoint[i].size()+1;
	}
	idx += gidx;
	return idx;
}
inline void FreePeak::treeIdx(const int bidx, int &tidx,int &gidx)const{
	int idx = bidx;
	tidx = 0, gidx = bidx;
	for (tidx = 0; tidx < m_tree.size(); ++tidx){
		idx -= m_tree[tidx]->region.size();
		if (idx < 0){			
			break;
		}
		gidx -= m_tree[tidx]->region.size();
	}	
}
inline int FreePeak::treeIdx(const int bidx)const{
	int idx = bidx,tidx = 0;
	for (tidx = 0; tidx < m_tree.size(); ++tidx){
		idx -= m_tree[tidx]->region.size();
		if (idx < 0)		break;		
	}
	return tidx;
}
inline int FreePeak::smallestBox(const set<int>& except){
	int sidx,tidx;
	double rat = 1;
	for (int i = 0; i < m_tree.size(); ++i){
		int idx = m_tree[i]->smallestBox();
		int bidx = boxIdx(i, idx);
		if (except.find(bidx) != except.end()) continue;
		double rs = m_tree[i]->region[idx].rat;
		if (rat > rs){
			tidx = i; sidx = idx;
		}
	}
	return boxIdx(tidx, sidx);
}
inline int FreePeak::largestBox(const set<int>& except){
	int sidx, tidx;
	double rat = 0.;
	for (int i = 0; i < m_tree.size(); ++i){
		int idx = m_tree[i]->largestBox();
		int bidx = boxIdx(i, idx);
		if (except.find(bidx) != except.end()) continue;
		double rs = m_tree[i]->region[idx].rat;
		if (rat < rs){
			tidx = i; sidx = idx;
		}
	}
	return boxIdx(tidx, sidx);
}
#endif
