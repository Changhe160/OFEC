#include "FSchwefel_2_13.h"


FSchwefel_2_13::FSchwefel_2_13(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
     mpp_a = new int * [m_numDim];
     mpp_b = new int * [m_numDim];
     mp_alpha = new double[m_numDim];
    for(int i=0;i<m_numDim; ++i) {
          mpp_a[i] = new int [m_numDim];
          mpp_b[i] = new int [m_numDim];
     }
     setSearchRange(-OFEC_PI, OFEC_PI);
     initialize();
}
FSchwefel_2_13::FSchwefel_2_13(const int rId, const int rDimNumber, string& rName):Problem(rId, rDimNumber, rName,1),\
	BenchmarkFunction(rId, rDimNumber, rName,1){
     mpp_a = new int * [m_numDim];
     mpp_b = new int * [m_numDim];
     mp_alpha = new double[m_numDim];
    for(int i=0;i<m_numDim; ++i) {
          mpp_a[i] = new int [m_numDim];
          mpp_b[i] = new int [m_numDim];
     }
     setSearchRange(-OFEC_PI, OFEC_PI);
     initialize();
}

FSchwefel_2_13::~FSchwefel_2_13(){
    //dtor
    for(int i=0;i<m_numDim; ++i) {
          delete[] mpp_a[i];
          delete[] mpp_b[i];
    }
    delete[] mpp_a;
    delete[] mpp_b;
    delete[] mp_alpha;
	mpp_a=0;
	mpp_b=0;
	mp_alpha=0;
}

void FSchwefel_2_13::initialize(){
    setOriginalGlobalOpt();
    loadData(); 
    setBias(-460);

	vector<double> v(m_numDim,0);
	copy(mp_alpha,mp_alpha+m_numDim,v.begin());
    setGlobalOpt(0,&v);

    setAccuracy(1.0e-2);
}

void FSchwefel_2_13::loadData()
{
     string sa;
	char astr[100];
	sprintf(astr,"%d",m_numDim);
	strcat(astr,"Dim.txt");
	sa=astr;
	sa.insert(0,m_name+"_a_");
	
	sa.insert(0,"Problem/FunctionOpt/Data/" );
	sa.insert(0,Global::g_arg[param_workingDir] );//probDataPath

	string sb;
	sprintf(astr,"%d",m_numDim);
	strcat(astr,"Dim.txt");
	sb=astr;
	sb.insert(0,m_name+"_b_");
	sb.insert(0,"Problem/FunctionOpt/Data/" );
	sb.insert(0,Global::g_arg[param_workingDir] );//probDataPath

	string salpha;
	sprintf(astr,"%d",m_numDim);
	strcat(astr,"Dim.txt");
	salpha=astr;
	salpha.insert(0,m_name+"_alpha_");
	
	salpha.insert(0,"Problem/FunctionOpt/Data/" );
	salpha.insert(0,Global::g_arg[param_workingDir] );//probDataPath

	ifstream in_a;
	in_a.open(sa.data());
	ifstream in_b;
	in_b.open(sb.data());
	ifstream in_alpha;
	in_alpha.open(salpha.data());
	if(in_a.fail()){
	     for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim; ++j) {
                    mpp_a[i][j] = int(-100.0 + Global::msp_global->mp_uniformPro->Next() * 200);
               }
	     }
          ofstream out(sa.c_str());
          for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim;j++)  {
                    out<<mpp_a[i][j]<<" ";
               }
          }
          out.close();
	}else{
	     for(int i=0;i<m_numDim; ++i) {
             for(int j=0;j<m_numDim;j++) {
                 in_a>>mpp_a[i][j];
               }
	     }
	}
	in_a.close();

	if(in_b.fail()){
	     for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim; ++j) {
                    mpp_b[i][j] = int(-100.0 + Global::msp_global->mp_uniformPro->Next() * 200);
               }
	     }
          ofstream out(sb.c_str());
          for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim;j++)  {
                    out<<mpp_b[i][j]<<" ";
               }
          }
          out.close();
	}else{
	     for(int i=0;i<m_numDim; ++i) {
             for(int j=0;j<m_numDim;j++) {
                 in_b>>mpp_b[i][j];
               }
	     }
	}
	in_b.close();

	if(in_alpha.fail()){
	     for(int i=0;i<m_numDim; ++i) {
               mp_alpha[i] = -OFEC_PI + Global::msp_global->mp_uniformPro->Next() * 2 * OFEC_PI;
	     }
          ofstream out(salpha.c_str());
          for(int i=0;i<m_numDim; ++i) {
                    out<<mp_alpha[i]<<" ";
          }
          out.close();
	}else{
	     for(int i=0;i<m_numDim; ++i) {
                 in_alpha>>mp_alpha[i];
	     }
	}
	in_alpha.close();
}



void FSchwefel_2_13::evaluate__(double const *x,vector<double>& obj){
     double result = 0;
     for(int i=0;i<m_numDim; ++i) {
          double A = 0;
          double B = 0;
          for(int j=0;j<m_numDim; ++j) {
               A += mpp_a[i][j] * sin(mp_alpha[j]) + mpp_b[i][j] * cos(mp_alpha[j]);
               B += mpp_a[i][j] * sin(x[j]) + mpp_b[i][j] * cos(x[j]);
          }
          result += pow((A - B), 2.0);
     }
     obj[0]= result + m_bias;
}
