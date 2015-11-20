#include "DTLZ1.h"

DTLZ1::DTLZ1(ParamMap &v) :Problem((v[param_proId]), (v[param_numObj] + v[param_interTest1] - 1), (v[param_proName]), v[param_numObj]), \
DTLZ((v[param_proId]), (v[param_numObj] + v[param_interTest1] - 1), (v[param_proName]), v[param_numObj])
{
	//default value m_k=5;
}

void DTLZ1::evaluate__(double const *x,vector<double>& obj)
{
//	if (x.size() != M_ + k_ - 1) return false; // #variables does not match
	
	double g = 0;
	for (int i = m_numObj-1; i < m_numDim; i += 1)
	{
		g += (x[i]-0.5)*(x[i]-0.5) - cos(20*OFEC_PI*(x[i]-0.5));
	}
	g = (m_numDim+1-m_numObj + g)*100;

	for (int m = 0; m < m_numObj; m += 1)
	{
		double product = 0.5*(1+g);
		int i = 0;
		for (; m_numObj >= 2+m && i <= m_numObj-2-m; i += 1)
			product *= x[i];
		if (m > 0)
			product *= (1 - x[i]);
		obj[m] = product;
	}
}



