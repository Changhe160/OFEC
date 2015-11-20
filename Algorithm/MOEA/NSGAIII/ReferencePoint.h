#ifndef REFERENCE_POINT__
#define REFERENCE_POINT__

/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 11 Jan 2015
// Last modified:

#include "../../../Utility/include.h"
// ----------------------------------------------------------------------------------
//		CReferencePoint
//
// Reference points play very important roles in NSGA-III. Individuals in the population
// are associated with reference points, and the survivors in the environmental selection
// are determined based on the niche count of the reference points.
//
// Check Algorithms 1-4 in the orignal paper for the usage of reference points.
// ----------------------------------------------------------------------------------

class ReferencePoint
{
public:
	explicit ReferencePoint(size_t s): mv_position_(s), m_member_size_(0) {}
	ReferencePoint(){}
	~ReferencePoint();

	const vector<double> & pos() const { return mv_position_; }
	vector<double> & pos() { return mv_position_; }

	size_t MemberSize() const { return m_member_size_; }
	bool HasPotentialMember() const { return !mv_potential_members_.empty(); }
	void clear();
	void AddMember();
	void AddPotentialMember(size_t member_ind, double distance);
	int FindClosestMember() const;
	int RandomMember() const;
	void RemovePotentialMember(size_t member_ind);

private:
	vector<double> mv_position_;

	// pair<indices of individuals in the population, distance>
	// note. only the data of individuals in the last considered front
	// will be stored.
	vector< pair<size_t, double> > mv_potential_members_; 
	size_t m_member_size_; 
};

// ----------------------------------------------------------------------------------
// GenerateReferencePoints():
//
// Given the number of objectives (M) and the number of divisions (p), generate the set of 
// reference points. Check Section IV-B and equation (3) in the original paper.

void GenerateReferencePoints(vector<ReferencePoint> *rps, size_t M, const vector<size_t> &p);
// ----------------------------------------------------------------------------------
// Associate():
//
// Associate individuals in the population with reference points.
// Check Algorithm 3 in the original paper.

void Associate(vector<ReferencePoint> *prps, const vector<vector<double> > &conv_obj, const vector<vector<int> > &fronts);
// ----------------------------------------------------------------------------------

#endif