#include "FShubert.h"


FShubert::FShubert(ParamMap &v):Problem((v[param_proId]), v[param_numDim],(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), v[param_numDim],(v[param_proName]),1){
	
	initialize();
}
FShubert::FShubert(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId,rDimNumber , rName,1){
	
	initialize();
}

FShubert::~FShubert(){
    //dtor
}
void FShubert::initialize(){
	vector<double> lower, upper;
	for(int i=0;i<m_numDim;i++){
		upper.push_back (10);
		lower.push_back (-10);
	}
	setSearchRange(lower,upper);
	setDisAccuracy(0.5);
	
	m_globalOpt.setNumOpts(m_numDim*pow(3,m_numDim));
	if(m_numDim==2){
		setAccuracy(1.e-3);
		m_globalOpt.setFlagLocTrue();		
		ifstream in;
		stringstream ss;
		ss<<Global::g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/"<<m_name<<"_Opt_"<<m_numDim<<"Dim.txt";
		in.open(ss.str().c_str());
		if(in.fail()){
			throw myException("cannot open data file@FShubert::initialize()");
		}
		for(int i=0;i<18;++i){
			double x0,x1;
			in>>x0>>x1;
			m_globalOpt[i].data().m_x[0]=x0; m_globalOpt[i].data().m_x[1]=x1; m_globalOpt[i].data().m_obj[0]=-186.7309;
		}
		in.close();
	}else if(m_numDim==3){
		setAccuracy(1.e-2);
		m_globalOpt.flagLoc()=false;
		m_globalOpt.flagGloObj()=true;
		m_globalOpt.setGloObj(vector<vector<double>>(m_numDim*pow(3,m_numDim),vector<double>(1,-2709.09)));
	}else if(m_numDim==4){
		setAccuracy(1.e-1);
		m_globalOpt.flagLoc()=false;
		m_globalOpt.flagGloObj()=true;
		m_globalOpt.setGloObj(vector<vector<double>>(m_numDim*pow(3,m_numDim),vector<double>(1,-39303.6))); 
	}else{
		m_globalOpt.flagLoc()=false;
		m_globalOpt.flagGloObj()=false;
	}
	m_originalGlobalOpt=m_globalOpt;
	addProTag(MMP);
}
void FShubert::evaluate__(double const *x,vector<double>& obj){
	double s=1;
	for(int j=0;j<m_numDim;++j){
		double a=0;
		for(int i=1;i<=5;i++){
			a+=i*cos((i+1)*x[j]+i);		
		}
		s*=a;
	}
	obj[0]= s+m_bias;

}