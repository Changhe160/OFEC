#include "ReferencePoint.h"
#include "../../../Global/global.h"
#include "MathAux.h"

#include <limits>
using namespace std;


ReferencePoint::~ReferencePoint()
{
	mv_potential_members_.clear();
	mv_position_.clear();
}

void ReferencePoint::clear()
{
	m_member_size_ = 0;
	mv_potential_members_.clear();
	
}
// ----------------------------------------------------------------------
void ReferencePoint::AddMember()
{
	m_member_size_ += 1;
}
// ----------------------------------------------------------------------
void ReferencePoint::AddPotentialMember(size_t member_ind, double distance)
{
	mv_potential_members_.push_back(move(make_pair(member_ind, distance)));
}
// ----------------------------------------------------------------------
int ReferencePoint::FindClosestMember() const
{
	double min_dist = numeric_limits<double>::max();
	int min_indv = -1;
	for (size_t i=0; i<mv_potential_members_.size(); i+=1)
	{
		if (mv_potential_members_[i].second < min_dist)
		{
			min_dist = mv_potential_members_[i].second;
			min_indv = mv_potential_members_[i].first;
		}
	}
	return min_indv;
}
// ----------------------------------------------------------------------
int ReferencePoint::RandomMember() const
{
	if (mv_potential_members_.size() > 0)
	{
		return mv_potential_members_[Global::msp_global->getRandInt(0,mv_potential_members_.size())].first;
	}
	else
	{
		return -1;
	}
}
// ----------------------------------------------------------------------
void ReferencePoint::RemovePotentialMember(size_t member_ind)
{
	for (size_t i=0; i<mv_potential_members_.size(); i+=1)
	{
		if (mv_potential_members_[i].first == member_ind)
		{
			mv_potential_members_.erase(mv_potential_members_.begin()+i);
			return;
		}
	}
}


// ----------------------------------------------------------------------
// Other utility functions
// ----------------------------------------------------------------------
void generate_recursive(vector<ReferencePoint> *rps, ReferencePoint *pt, size_t num_objs, 
						size_t left, size_t total, size_t element)
{
	if (element == num_objs-1)
	{
		pt->pos()[element] = static_cast<double>(left)/total;
		rps->push_back(*pt);
	}
	else
	{
		for (size_t i=0; i<=left; i+=1)
		{
			pt->pos()[element] = static_cast<double>(i)/total;
			generate_recursive(rps, pt, num_objs, left-i, total, element+1);
		}
	}
}
// ----------------------------------------------------------------------
void GenerateReferencePoints(vector<ReferencePoint> *rps, size_t M, const vector<size_t> &p)
{
	ReferencePoint pt(M);

	generate_recursive(rps, &pt, M, p[0], p[0], 0);

	if (p.size()>1) // two layers of reference points (Check Fig. 4 in NSGA-III paper)
	{
		vector<ReferencePoint> inside_rps;
		generate_recursive(&inside_rps, &pt, M, p[1], p[1], 0);

		double center = 1.0/M;

		for (size_t i=0; i<inside_rps.size(); i+=1)
		{
			for (size_t j=0; j<inside_rps[i].pos().size(); j+=1)
			{
				inside_rps[i].pos()[j] = (center+inside_rps[i].pos()[j])/2; // (k=num_divisions/M, k, k, ..., k) is the center point
			}
			rps->push_back(inside_rps[i]);
		}
	}
}
// ----------------------------------------------------------------------
void Associate(std::vector<ReferencePoint> *prps, const vector<vector<double> > &conv_obj, const vector<vector<int> > &fronts)
{
	vector<ReferencePoint> &rps = *prps;

	for (size_t t=0; t<fronts.size(); t+=1)
	{
		for (size_t i=0; i<fronts[t].size(); i+=1)
		{
			size_t min_rp = rps.size();
			double min_dist = numeric_limits<double>::max();
			for (size_t r=0; r<rps.size(); r+=1)
			{
				double d = MathAux::PerpendicularDistance(rps[r].pos(), conv_obj[ fronts[t][i] ]);
				if (d < min_dist)
				{
					min_dist = d;
					min_rp = r;
				}
			}

			if (t+1 != fronts.size()) // associating members in St/Fl (only counting)
			{
				rps[min_rp].AddMember();
			}
			else
			{
				rps[min_rp].AddPotentialMember(fronts[t][i], min_dist);
			}

		}// for - members in front
	}// for - fronts
}
// ----------------------------------------------------------------------