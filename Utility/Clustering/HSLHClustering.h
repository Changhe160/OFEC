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
// Last modified:
#ifndef HSLHCLUSTERING_H
#define HSLHCLUSTERING_H
// heuristic single-linkage hierarchical clustering
#include "Group.h"

#ifdef DEMON_OFEC
#include "../../buffer/buffer.h"
#endif

template<typename ED, typename TypeIndi=Individual< Solution<ED> > >
class Cluster{
protected:
	vector<Group<ED,TypeIndi>> m_group;
	int m_initialNumber;
	/*distance between all objects*/
	double **mpp_dis;
	double m_interDiversity,m_intraDiversity;
	/*distance between groups*/
	double **mpp_groupDis;
	double **mpp_normalGroupDis;
	// 17/02/2014
	int m_space;			//distance in solution space, objective space, or both
	bool m_normalFlag;		//flag of distance normalization
	int m_minGroupSize;		//minimum number of individuals of the largest group
protected:
	virtual void allocateMemory(const int rNum, const int rInitialNum);
	virtual void freeMemory();
	void computeDistMatrix();
	double computerGroupDist(const int idFrom,const int idTo);  //group ID
	bool isNearestGroup(int & g1, int & g2, const unsigned subSize, bool constrain=true);
	void deleteGroup(const int index, const int id=-1);
	void addGroup(Group<ED,TypeIndi> &g);
	void updateDiversity();
	void updateGroupDistance(const int id);
	double getGroupDistance(const int idxFrom, const int idxTo);
	//16/01/2014
	void nomalizeGroupDis();
public:
	Cluster( );
	Cluster(const TypeIndi *p,const int num);
	void initialize(const Cluster<ED,TypeIndi> &clst);

	template<typename Type>
	void initialize(Type &p,const int num);
	~Cluster(){
		freeMemory();
	}
	Cluster<ED,TypeIndi> &operator =(const Cluster<ED,TypeIndi> &clst);	
	void roughClustering(const int subSize);	
	//17/05/2013
	void adaptiveClustering();	
	//17/02/2014
	void setSpace(const int space){
		m_space=space;
	}
	void setNormalizationFlag(const bool flag){
		m_normalFlag=flag;
	}
	const vector<Group<ED,TypeIndi>> & getGroup(){
		return m_group;
	} 
	Group<ED,TypeIndi> & operator[](const int i){
		return m_group[i];
	}
	int getSize(){
		return m_group.size();
	}
	void clear_(){
		freeMemory();
	}
	int getTotalNumItems(){
		int num=0;
		for(auto &i:m_group) num+=i.getSize();
		return num;
	}
	void setMinGroupSize(int size){
		m_minGroupSize=size;
	}
	int getGourpSizeAbove(){
		int count=0;
		for(int i=0;i<m_group.size();i++){
			if(m_group[i].getSize()>=m_minGroupSize) ++count;
		}
		return count;
	}
};

template <typename ED,typename TypeIndi>
Cluster<ED,TypeIndi>::Cluster():m_initialNumber(0),m_group(),mpp_dis(0),mpp_groupDis(0),mpp_normalGroupDis(0),\
	m_space(0),m_normalFlag(false),m_interDiversity(0),m_intraDiversity(0){

}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::allocateMemory(const int rNum, const int rInitialNum){
	m_group.clear();
	//for(int i=0;i<rNum;i++) m_group.push_back(move(Group<TypeIndi>()));
	m_group.resize(rNum);

    mpp_dis=new double*[rInitialNum];
	for(int i=0;i<rInitialNum;i++)
		mpp_dis[i]=new double[rInitialNum];
	mpp_groupDis=new double*[rInitialNum];
	for(int i=0;i<rInitialNum;i++)
		mpp_groupDis[i]=new double[rInitialNum];

	mpp_normalGroupDis=new double*[rInitialNum];
	for(int i=0;i<rInitialNum;i++)
		mpp_normalGroupDis[i]=new double[rInitialNum];
}

template <typename ED,typename TypeIndi>
Cluster<ED,TypeIndi>::Cluster(const TypeIndi * p, const int num):m_initialNumber(num){
	allocateMemory(num,m_initialNumber);

	for(int i=0;i<num;i++)
		m_group[i].initialize(p[i],i);

	computeDistMatrix();
}
template <typename ED,typename TypeIndi> 
template<typename Type>
void Cluster<ED,TypeIndi>::initialize(Type &p,const int num){
	freeMemory();
	m_initialNumber=num;
	allocateMemory(num,m_initialNumber);

	for(unsigned i=0;i<m_group.size();i++)
		m_group[i].initialize(p[i],i);
	computeDistMatrix();
	updateDiversity();
}
template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::initialize(const Cluster<ED,TypeIndi> &clst){
	if(clst.m_initialNumber==0) return;

	m_initialNumber=clst.m_initialNumber;
	allocateMemory(clst.m_group.size(),m_initialNumber);

	for(int i=0;i<m_group.size();i++)
		m_group[i].initialize(clst.m_group[i]);

	for(int i=0;i<m_initialNumber;i++)
		for(int j=0;j<m_initialNumber;j++)
			mpp_dis[i][j]=clst.mpp_dis[i][j];

	for(int i=0;i<m_initialNumber;i++)
		for(int j=0;j<m_initialNumber;j++)
			mpp_groupDis[i][j]=clst.mpp_groupDis[i][j];

	for(int i=0;i<m_initialNumber;i++)
		for(int j=0;j<m_initialNumber;j++)
			mpp_normalGroupDis[i][j]=clst.mpp_normalGroupDis[i][j];
	m_space=clst.m_space;
	m_normalFlag=clst.m_normalFlag;
	m_interDiversity=clst.m_interDiversity;
	m_intraDiversity=clst.m_intraDiversity;
}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::computeDistMatrix(){
	// only called by constructor
	for(int i=0;i<m_initialNumber;i++){
		mpp_dis[i][i]=-1; //infinite large
		for(int j=0;j<i;j++){
			mpp_dis[i][j]=mpp_dis[j][i]=m_group[i].m_center.getDistance(m_group[j].m_center);
			mpp_groupDis[i][j]=mpp_groupDis[j][i]=mpp_dis[i][j];
		}
	}
	nomalizeGroupDis();
}
template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::roughClustering(const int subSize){
	while(1){
		decltype( m_group.size()) i=0;
		while(i<m_group.size()&&m_group[i].m_member.size()>1) i++;
		if(i==m_group.size()) break;
		int g1=0,g2=0;
		if(!isNearestGroup(g1,g2,subSize)) break;
		m_group[g2].merge(m_group[g1]);
		// corrected 18/02/2014
		int idg2=m_group[g2].m_ID;
		deleteGroup(g1);
		updateGroupDistance(idg2);
	}
	for(auto &i:m_group) i.calculateRadius();
}
template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::deleteGroup(const int index, const int ID){
	//delete either by index or ID
	if(m_group.size()<2){
		freeMemory();
		return;
	}
	//20/05/2013

	int idx=-1;
	if(index!=-1) idx=index;
	else	while(m_group[++idx].m_ID!=ID&&idx<getSize()); 

	m_group.erase(m_group.begin()+idx);	
}
template <typename ED,typename TypeIndi>
double Cluster<ED,TypeIndi>::computerGroupDist(const int idFrom,const int idTo){
	if(idFrom==idTo) return -1;

	int idxFrom=0,idxTo=0;
	for(unsigned i=0;i<m_group.size();i++){
		if(m_group[i].m_ID==idFrom) idxFrom=i;
		if(m_group[i].m_ID==idTo) idxTo=i;
	}
	return m_group[idxFrom].m_center.getDistance(m_group[idxTo].m_center);
}
template <typename ED,typename TypeIndi>
bool Cluster<ED,TypeIndi>::isNearestGroup(int & g1, int & g2,const unsigned subSize, bool constrain){
	
	if(m_group.size()<=1) {
	 g1=0;
	 g2=0;
	 return false;
	}
	bool flag_fail=true;
	double Min_dis=0,dist;

	Min_dis=getGroupDistance(0,1);

	for(unsigned i=0;i<m_group.size();i++){
		// can't merge two mp_groups whose m_number are both greater than g_subSize
		for(unsigned j=0;j<m_group.size();j++){
				if(j==i) continue;
				if(constrain&&m_group[i].m_member.size()+m_group[j].m_member.size()>subSize) continue;
				//dist=computerGroupDist(i,j);
				dist=getGroupDistance(i,j);
				if(Min_dis>dist){
					Min_dis=dist;
					g1=i;
					g2=j;
					flag_fail=false;
				}
		}
	}

	return !flag_fail;

}


template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::freeMemory(){
	if(m_group.size()>0){
		m_group.clear();
		for(int i=0;i<m_initialNumber;i++){
			delete [] mpp_dis[i];
			mpp_dis[i]=0;
			delete [] mpp_groupDis[i];
			mpp_groupDis[i]=0;
			
			delete [] mpp_normalGroupDis[i];
			mpp_normalGroupDis[i]=0;
		}
		delete[] mpp_dis;
		delete [] mpp_groupDis;
		delete [] mpp_normalGroupDis;
		mpp_dis=0;

		m_initialNumber=0;
		mpp_groupDis=0;
		mpp_normalGroupDis=0;
	}
}
template <typename ED,typename TypeIndi>
Cluster<ED,TypeIndi> & Cluster<ED,TypeIndi>::operator =(const Cluster<ED,TypeIndi> &clst){
	if(clst.m_initialNumber==0) return *this;
	
	freeMemory();
	initialize(clst);

	return *this;
}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::adaptiveClustering(){
	if(m_group.size()<m_minGroupSize) {
		throw myException("The number of items is too few@Cluster<ED,TypeIndi>::adaptiveClustering()");
	}

	while(1){
		int g1=0,g2=0;
		int num=0;
		for(unsigned i=0;i<m_group.size();i++) {
			if(m_group[i].m_member.size()>1) num++;
		}
		if(num==m_group.size()) break;
		if(!isNearestGroup(g1,g2,0,false)) break;

		m_group[g2].merge(m_group[g1]);
		int idg2=m_group[g2].m_ID;
		deleteGroup(g1);
		updateGroupDistance(idg2);		
		updateDiversity();

		if(m_interDiversity<=m_intraDiversity&&getGourpSizeAbove()>0){
			break;
		}

	}

}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::updateDiversity(){
	m_interDiversity=m_intraDiversity=0;
	//calculate inter-diversity
	for(unsigned i=0;i<m_group.size();i++){
		for(unsigned j=0;j<i;j++){
			m_interDiversity+=mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID];
		}
	}
	// calulate intra-diversity
	for(unsigned i=0;i<m_group.size();i++){
		double div=0;
		for(unsigned j=0;j<m_group[i].m_member.size();j++){
			for(unsigned k=0;k<j;k++){
				div+=mpp_dis[m_group[i].m_member[j]->getId()-1][m_group[i].m_member[k]->getId()-1];
			}
		}
		m_intraDiversity+=div;
	}
	
	//Global::g_interDiversity=m_interDiversity;
	//Global::g_intraDiversity=m_intraDiversity;
}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::updateGroupDistance(const int id){
	int objID=-1;
	for(unsigned i=0;i<m_group.size();i++) {
		objID=m_group[i].m_ID;
		mpp_groupDis[objID][id]=mpp_groupDis[id][objID]=computerGroupDist(objID,id);
	}	
	nomalizeGroupDis();
}

template <typename ED,typename TypeIndi>
double Cluster<ED,TypeIndi>::getGroupDistance(const int idxFrom, const int idxTo){
	
	if(!m_normalFlag)	return mpp_groupDis[m_group[idxFrom].m_ID][m_group[idxTo].m_ID];
	else return mpp_normalGroupDis[m_group[idxFrom].m_ID][m_group[idxTo].m_ID];
	
}



template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::addGroup(Group<ED,TypeIndi> &gr){
	m_group.push_back(gr);
}

template <typename ED,typename TypeIndi>
void Cluster<ED,TypeIndi>::nomalizeGroupDis(){
	switch(m_space){
		case 0:{// solution space
			double maxDis=0;
			for(unsigned i=0;i<m_group.size();i++){
				for(unsigned j=0;j<i;j++){
					if(maxDis<mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID]) maxDis=mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID];
				}
			}
			for(unsigned i=0;i<m_group.size();i++){
				mpp_normalGroupDis[m_group[i].m_ID][m_group[i].m_ID]=-1;
				for(unsigned j=0;j<i;j++){
					mpp_normalGroupDis[m_group[i].m_ID][m_group[j].m_ID]=mpp_normalGroupDis[m_group[j].m_ID][m_group[i].m_ID]=mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID]/maxDis;
				}
			}
			break;
			}
		case 1:// objective space
			break;
		case 2:{//solution & objective space
			double maxDis=0;
			for(unsigned i=0;i<m_group.size();i++){
				for(unsigned j=0;j<i;j++){
					if(maxDis<mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID]) maxDis=mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID];
				}
			}
			for(unsigned i=0;i<m_group.size();i++){
				for(unsigned j=0;j<m_group.size();j++){
					if(i==j) mpp_normalGroupDis[m_group[i].m_ID][m_group[j].m_ID]=-1;
					else	mpp_normalGroupDis[m_group[i].m_ID][m_group[j].m_ID]=mpp_groupDis[m_group[i].m_ID][m_group[j].m_ID]/maxDis;
				}
			}
			double maxFit,minFit;
			double *normalFit=new double [m_group.size()];
			maxFit=minFit=m_group[0].m_center.obj(0);
			for(unsigned i=1;i<m_group.size();i++){
				normalFit[i]=m_group[i].m_center.obj(0);
				if(maxFit<normalFit[i]) maxFit=normalFit[i];
				if(minFit>normalFit[i]) minFit=normalFit[i];
			}
			double disFit=maxFit-minFit;
			for(unsigned i=0;i<m_group.size();i++){
				normalFit[i]=(normalFit[i]-minFit)/ disFit;
			}
			
			//2 times normalized distance in solution space + normalized distance in objective space
			for(unsigned i=0;i<m_group.size();i++){
				for(unsigned j=0;j<m_group.size();j++){
					if(i!=j)
					mpp_normalGroupDis[m_group[i].m_ID][m_group[j].m_ID]=3*mpp_normalGroupDis[m_group[i].m_ID][m_group[j].m_ID]+fabs(normalFit[i]-normalFit[j]);
				}
			}

			delete [] normalFit;
			normalFit=0;
			break;
		}	
	}
}

#endif // CLUSTERING_H
