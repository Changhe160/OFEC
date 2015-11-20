#include "DTLZ4.h"

DTLZ4::DTLZ4(ParamMap &v):Problem((v[param_proId]), (v[param_numObj]+v[param_interTest1]-1),(v[param_proName]),v[param_numObj]),\
DTLZ((v[param_proId]), (v[param_numObj] + v[param_interTest1] - 1), (v[param_proName]), v[param_numObj])
{
	//default value m_k=10;
}

void DTLZ4::evaluate__(double const *x,vector<double>& obj)
{
	double g = 0;
	for (size_t i = m_numObj-1; i < m_numDim; i += 1)
		g += pow((x[i]-0.5),2);

	for (size_t m = 0; m < m_numObj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=m_numObj-2; i+=1)
			product *= cos( pow(x[i],100)*OFEC_PI/2);
		if (m > 0)
			product *= sin( pow(x[i],100)*OFEC_PI/2);		
		obj[m] = product;
	}
}



