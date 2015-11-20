#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "../Utility/include.h"
#include "../Utility/TypeList/Typelist.h"

struct DomainReal{
		double m_upper;					//Dimension boudary in the landscape
		double m_lower;					//**************************************************//
		bool m_flag;                        // flag of rigid boundary
};
class BoundaryCont{						//***************************************************//

protected:
	unsigned  m_numDim;
	double m_domainSize;
	vector<DomainReal> m_domain;
protected:
	virtual void calculateDomainSize();
public:
	BoundaryCont(int numDim=1);
	void setBoundary(const vector<double>& rLower, const vector<double> &rUpper);
	void setBoundary(double rLower, double rUpper, int i=-1);
	void setBoundaryFlag(bool flag);
	void setBoundaryFlag(const vector<bool> &flag);
	bool getBoundaryFlag(int i);
	BoundaryCont & operator=(const BoundaryCont &rBound);
	double getDomainSize();
	DomainReal &operator[](const int i);
	const DomainReal &operator[](const int i)const;
	void getSearchRange(double &l,double &u, int i);
	DomainReal & getDomain(int i);
	double getMinRange();
	void resizeDim(const int num);
};
inline void BoundaryCont::getSearchRange(double &l, double &u, int i){
	l = m_domain[i].m_lower;
	u = m_domain[i].m_upper;
}

struct BoundaryTSP{
	int m_upper;					
	int m_lower;
	void setBoundary(int l, int u){
		m_lower=l; m_upper=u;
	}
};

enum BoundaryType{Bound_Real=0,Bound_TSP};

typedef LOKI_TYPELIST_2(BoundaryCont,BoundaryTSP) BoundaryTypeList;

typedef Loki::TL::TypeAt<BoundaryTypeList,Bound_Real>::Result Boundary;

#endif