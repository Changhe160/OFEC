#ifndef PROG_RUN_H
#define PROG_RUN_H
#include "../Utility/TypeVar/typeVar.h"

struct Run
{
	int m_runId;
	Run(int runId,ParamMap &v);
	~Run();
	ReturnFlag go();
	void test();
};

#endif