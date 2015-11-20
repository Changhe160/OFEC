#include "ZDT6.h"

ZDT6::ZDT6(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),2),ZDT(v)
{
  generateAdLoadPF();
}

void ZDT6::evaluate__(double const *x,vector<double>& obj)
{
	double g = 0;
	for(int n=1;n<m_numDim;n++)
		g=g+x[n];
	g=pow(g/(m_numDim-1),0.25);
	g = 1 + 9*g;
	obj[0] = 1-exp(-4*x[0])*pow(sin(6*OFEC_PI*x[0]),6);
	obj[1] = g*(1 - pow(obj[0]/g,2));
}





