#include "OneMax.h"

void OneMax::setObjSet() {
	m_os.clear();
	if (!m_globalOpt.flagGloObj()) return;
	int num = m_globalOpt.getNumOpt();
	for (int i = 0; i < num; ++i) {
		m_os.push_back(&m_globalOpt[i].data().m_obj);
	}
}
OneMax::OneMax(ParamMap& v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1) {
	setOptType(MAX_OPT);
	addProTag(ONEMAX);
	m_globalOpt.setNumOpts(1);
	for (int i = 0; i < m_numDim; i++)	m_globalOpt[0].data().m_x[i] = 1;
	m_globalOpt[0].data().m_obj[0] = m_numDim;
}
OneMax::OneMax(const int rId, const int rDimNumber, string rName, const int numObj) :Problem(rId, rDimNumber, rName, numObj) {

	setOptType(MAX_OPT);
	addProTag(ONEMAX);
	m_globalOpt.setNumOpts(1);
	for (int i = 0; i < m_numDim; i++)	m_globalOpt[0].data().m_x[i] = 1;
	m_globalOpt[0].data().m_obj[0] = m_numDim;

}

ReturnFlag OneMax::evaluate_(VirtualEncoding &s_, bool rFlag, ProgramMode mode, bool flag2) {
	CodeVInt &s = dynamic_cast< CodeVInt&>(s_);

	for (int i = 0; i<m_numObj; i++)
		s.m_obj[i] = 0;

	for (int n = 0; n<m_numObj; n++)
	{
		for (size_t i = 0; i<m_numDim; i++)
		{
			s.m_obj[n] +=s.m_x[i];
		}
	}
	if (flag2) {
		if (rFlag)	m_evals++;

		if (Global::msp_global->mp_algorithm.get() != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()) return Return_Terminate;
		return Return_Normal;
	}
	return Return_Normal;
}

bool OneMax::isValid(const VirtualEncoding &s_) {
	const CodeVInt &s = dynamic_cast< const CodeVInt&>(s_);
	for (size_t i = 0; i<m_numDim; i++)
	{
		if (s.m_x[i] != 0 && s.m_x[i] != 1) return false;
	}
	return true;
}

void OneMax::initializeSolution(VirtualEncoding &result, const int idx, const int maxId) {
	CodeVInt &s = dynamic_cast<  CodeVInt&>(result);
	for (size_t i = 0; i<m_numDim; i++)
	{
		if (Global::msp_global->mp_uniformPro->Next() < 0.5) s.m_x[i] = 0;
		else s.m_x[i] = 1;
	}
}

double  OneMax::getDistance(const VirtualEncoding &s1, const VirtualEncoding &s2, DistanceMode mode) {
	const CodeVInt& ss1 = dynamic_cast<const CodeVInt&>(s1);
	const CodeVInt& ss2 = dynamic_cast<const CodeVInt&>(s2);
	int dis = 0;
	for (size_t i = 0; i<m_numDim; i++)
	{
		if (ss1.m_x[i]!=ss2.m_x[i]) dis++;
	}
	return dis;
}

bool OneMax::isSame(const VirtualEncoding &s1, const VirtualEncoding &s2) {
	const CodeVInt& ss1 = dynamic_cast<const CodeVInt&>(s1);
	const CodeVInt& ss2 = dynamic_cast<const CodeVInt&>(s2);
	for (size_t i = 0; i<m_numDim; i++)
	{
		if (ss1.m_x[i] != ss2.m_x[i]) return false;
	}
	return true;
}