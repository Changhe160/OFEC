#include "F2.h"

F2::F2(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]),2),\
	F_Base((v[param_proId]), (v[param_numDim]),(v[param_proName]),2)
{
	vector<double>l(m_numDim,-1), u(m_numDim,1);
	l[0] = 0;
	setSearchRange(l,u);
	m_dtype=1;
	m_ptype=21;
	m_ltype=22;
	LoadPF();
}

void F2::evaluate__(double const *x,vector<double>& obj)
{
	F_Base::evaluate__(x,obj);
}
