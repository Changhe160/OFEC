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
/*
Yuhui Shi; Eberhart, R., "A modified particle swarm optimizer," Evolutionary Computation 
Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE 
International Conference on , vol., no., pp.69,73, 4-9 May 1998
*/
#ifndef PARTICLE_H
#define PARTICLE_H

#include "../Individual.h"
#include "../../Problem/ContinuousProblem.h"

template <typename, typename>	class Swarm;
//inertia weight model
class Particle:public Individual<CodeVReal>{
	template <typename,typename>	friend class Swarm;
	struct VMax{
		double m_min,m_max;
		VMax():m_min(0),m_max(0){}
	};
protected:
    Solution<CodeVReal> m_pbest;
	vector<double> m_vel;
	vector<VMax> m_vMax;
public:
    Particle();
	virtual ~Particle();
	Particle( const Particle& other);
	Particle( const Solution<CodeVReal>& other);
	Particle& operator=(const Particle & other);

    void initializeVmax();

    double getVelocity();
    void initializeVelocity();

    ReturnFlag NormalMutation(double *avg_v);
    virtual ReturnFlag move( const Solution<CodeVReal> & lbest, const Solution<CodeVReal> &gbest,double w, double c1, double c2);
	ReturnFlag moveBound( const Solution<CodeVReal> & lbest,  const Solution<CodeVReal> &gbest,double w, double c1, double c2);
 
    virtual ReturnFlag initialize(bool mode=true);
    virtual ReturnFlag initialize( int idex, int id,bool mode=true);
    virtual void initialize(const Solution<CodeVReal> &p, int idex, int id);
    virtual ReturnFlag initialize(const Solution<CodeVReal> &p,double radius, int idex, int id,bool mode=true);
    virtual ReturnFlag initialize(int rIdx,int rID, int rPopsize,bool mode=true);
    virtual ReturnFlag initialize(Solution<CodeVReal> *w, int mode=1,bool mode2=true);
    virtual void printToFile(ofstream & out);
    Solution<CodeVReal> & representative();
	const Solution<CodeVReal> & representative() const;
    void increaseDimension();
    void decreaseDimension();
    void updateMemory();
	void setVmax(double rMin, double rMax);
	void setVmax(double *rMin, double *rMax);
	void printToScreen();
	vector<double> getVel() const { return m_vel; }
	bool isSame(const Particle &p)const;
};

#endif // PARTICLE_H
