#include "Termination.h"
#include "../Global/global.h"

#ifdef OFEC_DEMON
extern bool g_algTermination;
#endif

bool Termination::ifTerminating() {
#ifdef OFEC_CONSOLE
	if (Global::msp_global->mp_problem != nullptr) {
		if (!Global::msp_global->mp_problem->isProTag(DOP)&&Global::msp_global->mp_problem->isGlobalOptFound())
			return true;
	}
#endif

#ifdef OFEC_DEMON
	if (g_algTermination) return true;
#endif

	return false;
};

bool TermMaxFes::ifTerminating() {

#if defined OFEC_DEMON
	return Termination::ifTerminating();
#else
	if (Termination::ifTerminating()) return true;
#endif

	if (Global::msp_global->mp_problem->getEvaluations() >= m_maxFes) return true;

	return false;
}


bool TermMean::ifTerminating(double value) {
#if defined OFEC_DEMON
	return Termination::ifTerminating();
#else
	if (Termination::ifTerminating()) return true;
#endif

	if (m_sucIter >= m_maxSucIter) return true;
	return false;
}

void TermMean::countSucIter(double value) {
	m_curMean = value;
	if (m_curMean < 1) m_epsilon = 1.e-5;

	if (Global::msp_global->mp_problem->getOptType() == MIN_OPT) {
		if (m_preMean - m_curMean>m_epsilon) {
			m_sucIter = 0;
			m_preMean = m_curMean;
		}
		else	m_sucIter++;
	}
	else {
		if (m_curMean - m_preMean >m_epsilon) {
			m_sucIter = 0;
			m_preMean = m_curMean;
		}
		else	m_sucIter++;
	}
}