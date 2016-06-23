#ifndef CMAES_H
#define CMAES_H

#include <string>
using namespace std;

#include "../../Algorithm.h"
#include "../../../Global/solution.h"

class CMAES:public Algorithm
{
     public:
        CMAES(ParamMap &v);
        ~CMAES();
		ReturnFlag run_();
protected:
	void copy(double*x, CodeVReal &vx);
	double fitCompute(Solution<CodeVReal> &s);
     private:
          string signalsFilePathName;
          string initialsFilePathName;

		  vector<Solution<CodeVReal>> m_pop;
};

#endif // CMAES_H
