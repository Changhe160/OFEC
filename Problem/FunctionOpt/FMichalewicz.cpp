#include "FMichalewicz.h"

FMichalewicz::FMichalewicz(ParamMap &v):Problem((v[param_proId]), 2,(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), 2,(v[param_proName]),1),m_m(20){
	v[param_numDim]=2;

	initialize();
}
FMichalewicz::FMichalewicz(const int rId, const int rDimNumber, string& rName):Problem(rId, 2, rName,1),\
	BenchmarkFunction(rId, 2, rName,1),m_m(20){

	initialize();
}

FMichalewicz::~FMichalewicz(){
    //dtor
}
void FMichalewicz::initialize(){
	vector<double> lower, upper;
	upper.push_back (OFEC_PI);
	upper.push_back (OFEC_PI);
	lower.push_back (0);
	lower.push_back (0);
	setSearchRange(lower,upper);

	setAccuracy(1.e-3);
	setDisAccuracy(0.2);
	setOptType(MAX_OPT);
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(2); //1 gopt + 1 lopt
	CodeVReal x(m_numDim,1);
	x.m_x[0]=2.20291; x.m_x[1]=1.5708; x.m_obj[0]=1.8013;
	m_globalOpt[0].data()=x;
	x.m_x[0]=2.20291; x.m_x[1]=2.71157; x.m_obj[0]=1.21406;
	m_globalOpt[1].data()=x;

	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
	setObjSet();
}
void FMichalewicz::evaluate__(double const *x,vector<double>& obj){
	double s=0;
	for(int i=0;i<m_numDim;++i){
		s+=sin(x[i])*pow(sin((i+1)*x[i]*x[i]/OFEC_PI),m_m);
	}
	obj[0]= s+m_bias;
}