#ifndef FSCLUSTERING_H
#define FSCLUSTERING_H

// Clustering by fast search and find of density peaks
// Science 27 June 2014: 
// Vol. 344 no. 6191 pp. 1492-1496 
// DOI: 10.1126/science.1242072
// http://www.sciencemag.org/content/344/6191/1492.full// 
//
// Code Author: Eric Yuan
// Blog: http://eric-yuan.me
// You are FREE to use the following code for ANY purpose.
//
// Have fun with it
#include "../include.h"

class FSCluster
{
public:
	
	struct Result{
		struct Infor{
			int clstNo;	//cluster No. 
			bool isCenter,isCore;
		};
		vector<Infor> info;
		int numClst;
	};
	enum LocalDensity {Gaussian_kernel=0, CutOff_kernel=1}; 

	FSCluster(const vector<vector<double> > &dis);
	void setMinRhondMinDelta(double,double);
	void setDisMatrix(const vector<vector<double> > &);
	void clustering(FSCluster::Result &result);
private:
	void getdc();
	void getLocalDensity(LocalDensity LocalDensityVersion=Gaussian_kernel);
	void getDistanceToHigherDensity();
	void getClusterNum();
	void assign();
	void getHalo(FSCluster::Result &result);
	void autoSetMinRhondMinDelta(double val);
private:
	int m_numSamples;
	vector<double> mv_rho;
	vector<double> mv_delta;
	vector<int> mv_cl;
	vector<int> mv_nneigh;
	vector<int> mv_icl;
	vector<int> mv_halo;
	double m_dc,m_rhomin,m_deltamin;
	int m_NCLUST;
	vector<vector<double> > mv_dis;
	
};


#endif