#include "FSchwefel_2_6.h"

FSchwefel_2_6::FSchwefel_2_6(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),\
	BenchmarkFunction((v[param_proId]), (v[param_numDim]),(v[param_proName]),1){
     mpp_a = new int * [m_numDim];
    for(int i=0;i<m_numDim; ++i) {
          mpp_a[i] = new int [m_numDim];
     }
     mp_b = new double [m_numDim];
    setSearchRange(-100, 100);
    initialize();
}
FSchwefel_2_6::FSchwefel_2_6(const int rId, const int rDim, string& rName):Problem(rId, rDim, rName,1),\
	BenchmarkFunction(rId, rDim, rName,1){
     mpp_a = new int * [m_numDim];
    for(int i=0;i<m_numDim; ++i) {
          mpp_a[i] = new int [m_numDim];
     }
     mp_b = new double [m_numDim];
    setSearchRange(-100, 100);
    initialize();
}

FSchwefel_2_6::~FSchwefel_2_6(){
    //dtor
    for(int i=0;i<m_numDim; ++i) {
          delete[] mpp_a[i];
    }
    delete[] mpp_a;
    delete[] mp_b;
	mpp_a=0;
	mp_b=0;
}

void FSchwefel_2_6::initialize(){
     setOriginalGlobalOpt();
	 if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_2_6_Bound_CEC05")){
	 
	 }else{
		throw myException("Error: please check the problem ID@FSchwefel_2_6::initialize");
	 }

 
    loadData();
    setBias(-310);
    setGlobalOpt();
    setAccuracy(1.0e-2);
}

void FSchwefel_2_6::loadData()
{
    string sa;
	stringstream ss;
	ss<<m_numDim<<"Dim.txt";
	sa=ss.str();
	sa.insert(0,m_name+"_a_");
	sa.insert(0,"Problem/FunctionOpt/Data/" );
	sa.insert(0,Global::g_arg[param_workingDir] );//probDataPath

	ifstream in_a;
	in_a.open(sa.data());
	if(in_a.fail()){
	     for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim; ++j) {
                    mpp_a[i][j] = int(-500.0 + Global::msp_global->mp_uniformPro->Next() * 1000);
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


	string so;
	so=ss.str();
	so.insert(0,m_name+"_Opt_");

	so.insert(0,"Problem/FunctionOpt/Data/" );
	so.insert(0,Global::g_arg[param_workingDir] );//probDataPath

	ifstream in;
	in.open(so.data());
	if(in.fail()){
            for(int j=0;j<m_numDim;j++){
                double rl,ru,range;
				ru=m_searchRange[j].m_upper-(m_originalGlobalOpt[0].data().m_x[j]);
				rl=(m_originalGlobalOpt[0].data().m_x[j])-m_searchRange[j].m_lower;
                range=rl<ru?rl:ru;
                mp_translation[j]=(Global::msp_global->mp_uniformPro->Next()-0.5)*2*range;
            }
            for(int i=0;i<m_numDim; ++i) {
                    if(i < m_numDim / 4) mp_translation[i] = -100;
                    else if(i >= m_numDim * 3 / 4 - 1) mp_translation[i] = 100;
               }

          ofstream out(so.c_str());
          for(int j=0;j<m_numDim;j++)        out<<mp_translation[j]<<" ";
		out.close();
	}
	else{

        for(int j=0;j<m_numDim;j++) {
            in>>mp_translation[j];

		}
	}
	in.close();

	for(int i=0;i<m_numDim; ++i) {
          mp_b[i] = 0;
          for(int j=0;j<m_numDim; ++j) {
               mp_b[i] += mpp_a[i][j] * mp_translation[j];
          }
	}
}


void FSchwefel_2_6::evaluate__(double const *x,vector<double>& obj){

     double fit = 0;
     double * tempVector = new double[m_numDim];
	 if(IS_PROBLEM_NAME(m_id,"FUN_Schwefel_2_6_Bound_CEC05")) {
          for(int i=0;i<m_numDim; ++i) {
               for(int j=0;j<m_numDim; ++j) {
                    tempVector[j] = 0;
                    for(int k=0; k<m_numDim; ++k) {
                         tempVector[j] += mpp_a[j][k] * x[k];
                    }
               }
               for(int j=0;j<m_numDim; ++j) {
                    tempVector[j] -= mp_b[j];
               }
               double temp = 0;
               for(int  j=0; j<m_numDim; ++j) {
                    temp += pow(tempVector[j], 2.0);
               }
               temp = sqrt(temp);
               if(fit < temp) fit = temp;
          }
     }
     delete[] tempVector;
	 tempVector=0;
    obj[0]= fit+m_bias;
}
