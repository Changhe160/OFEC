#ifndef SYS_INITIALIZATION_H
#define SYS_INITIALIZATION_H
#include "../Utility/include.h"
#include "../Utility/TypeVar/typeVar.h"
struct Info
{
	string Name;
	int ClassIndex;
	set<string> Type;
};

void LoadData(vector<Info> & proInfo,vector<Info> & algInfo, vector<Info> &termInfo);
void registerProNdAlg();
void setAddtParameters();
void setGlobalParameters(int argn,char *argv[]);
void printGArg();
void setAlgParameters();

template<typename Base, typename Derived>
Base * createObject(ParamMap &v) {
	return new Derived(v);
}

#endif
