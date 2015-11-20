#include "F1.h"

F1::F1(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]),2),\
	F_Base((v[param_proId]), (v[param_numDim]),(v[param_proName]),2) 
{
	m_dtype=1;
	m_ptype=21;
	m_ltype=21;
	LoadPF();
}

void F1::evaluate__(double const *x,vector<double>& obj)
{
	F_Base::evaluate__(x,obj);
}
