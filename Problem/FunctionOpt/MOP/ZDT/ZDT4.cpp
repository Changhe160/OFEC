#include "ZDT4.h"

ZDT4::ZDT4(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),2),ZDT(v)
{
  generateAdLoadPF();
}

void ZDT4::evaluate__(double const *x,vector<double>& obj)
{
	double g = 0;
	for(int n=1;n<m_numDim;n++)
		g=g+(pow(x[n],2)-10*cos(4*OFEC_PI*x[n]));
	g = 1 + 10*(m_numDim-1)+g;
	obj[0] = x[0];
	obj[1] = g*(1 - sqrt(x[0]/g));
}





