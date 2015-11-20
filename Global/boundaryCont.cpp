#include "boundary.h"
#include "../Utility/myexcept.h"
#include "../Utility/definition.h"

void BoundaryCont::calculateDomainSize(){		
	m_domainSize=0.;
	for_each(m_domain.begin(),m_domain.end(),[this](DomainReal &d){	m_domainSize=m_domainSize+(d.m_upper-d.m_lower)*(d.m_upper-d.m_lower);});
	m_domainSize=sqrt(m_domainSize/m_numDim);	
}

BoundaryCont::BoundaryCont(int numDim):m_numDim(numDim),m_domain(numDim){
	for_each(m_domain.begin(),m_domain.end(), [](DomainReal &d){d.m_flag=true;});
}
void BoundaryCont::setBoundary(double rLower, double rUpper, int i){
	if(i>=(signed)m_domain.size()) throw myException("Out of range @ BoundaryCont::setBoundary(const TypeVar& rLower, const TypeVar &rUpper, int i)");
	
	if(i<0){
		for(auto &d:m_domain){ d.m_lower=rLower; d.m_upper=rUpper;}
	}else{
		m_domain[i].m_lower=rLower;
		m_domain[i].m_upper=rUpper;
	}
	calculateDomainSize();
}

void BoundaryCont::setBoundary(const vector<double>& rLower, const vector<double> &rUpper){
	for(decltype(m_domain.size()) i=0;i<m_domain.size();++i){
		m_domain[i].m_lower=rLower[i];
		m_domain[i].m_upper=rUpper[i];
	}
	calculateDomainSize();
}
void BoundaryCont::setBoundaryFlag(bool flag){ 
	for_each(m_domain.begin(),m_domain.end(), [&](DomainReal &d){d.m_flag=flag;});		
}
void BoundaryCont::setBoundaryFlag(const vector<bool> &flag){ 
	for(unsigned i=0;i<m_domain.size();i++){
		m_domain[i].m_flag=flag[i];
	}	
}
BoundaryCont & BoundaryCont::operator=(const BoundaryCont &rBound){
	if(this==&rBound) return *this;
	if(rBound.m_numDim!=m_numDim) throw myException("the number of dimensions must be the same@BoundaryCont::operator=");

	copy(rBound.m_domain.begin(),rBound.m_domain.end(), m_domain.begin());
	m_domainSize=rBound.m_domainSize;

	return *this;
};

double BoundaryCont::getDomainSize(){ 
		return m_domainSize; 
}
DomainReal & BoundaryCont::operator[](const int i){
	return m_domain[i];
}
const DomainReal &BoundaryCont::operator[](const int i)const{
	return m_domain[i];
}

bool BoundaryCont::getBoundaryFlag(int i){
	return m_domain[i].m_flag;
}

DomainReal & BoundaryCont::getDomain(int i){
	return m_domain[i];
}
double BoundaryCont::getMinRange(){
	double val = DBL_MAX;
	for (auto &i : m_domain){
		if (val > i.m_upper - i.m_lower)val = i.m_upper - i.m_lower;
	}
	return val;
}

void BoundaryCont::resizeDim(const int size){
	m_numDim = size;
	m_domain.resize(size);
	calculateDomainSize();
}