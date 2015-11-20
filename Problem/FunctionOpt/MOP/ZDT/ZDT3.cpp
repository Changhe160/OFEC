#include "ZDT3.h"

ZDT3::ZDT3(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),2),ZDT(v)
{
  generateAdLoadPF();
}

void ZDT3::evaluate__(double const *x,vector<double>& obj)
{
	double g = 0;
	for(int n=1;n<m_numDim;n++)
		g=g+x[n];
	g = 1 + 9*g/(m_numDim-1);
	obj[0] = x[0];
	obj[1] = g*(1 - sqrt(x[0]/g)-x[0]*sin(10*x[0]*OFEC_PI)/g);
}





