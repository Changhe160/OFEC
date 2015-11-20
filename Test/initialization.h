#ifndef SYS_INITIALIZATION_H
#define SYS_INITIALIZATION_H
#include "../Utility/include.h"
struct Info
{
	string Name;
	int ClassIndex;
	set<string> Type;
};

void LoadData(vector<Info> & proInfo,vector<Info> & algInfo);
void registerProNdAlg();
void setAddtParameters();
void initializeMeasueMulPop();
void setGlobalParameters(int argn,char *argv[]);
void printGArg();
void setAlgParameters();
#endif
