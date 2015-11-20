#include "F6.h"

F6::F6(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]),3),\
	F_Base((v[param_proId]), (v[param_numDim]),(v[param_proName]),3) 
{
	vector<double>l(m_numDim, -2), u(m_numDim, 2);
	l[0] = 0; u[0] = 1;
	l[1] = 0; u[1] = 1;
	setSearchRange(l, u);
	m_dtype=1;
	m_ptype=31;
	m_ltype=32;
	LoadPF();
}

void F6::evaluate__(double const *x,vector<double>& obj)
{
	F_Base::evaluate__(x,obj);
}