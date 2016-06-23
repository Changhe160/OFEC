#include "mMultiObj.h"
//#include "../Problem/FunctionOpt/MOP/FreePeak_M_OnePeak/FFreePeak_M_OnePeak.h"
unique_ptr<mMultiObj> mMultiObj::msp_perf(nullptr);

mMultiObj::mMultiObj(Global * glob,ParamMap &v)
{
	int maxNumrun=MAX_NUM_RUN;
	int numDim=glob->mp_problem->getNumDim();
	int numObj=glob->mp_problem->getNumObj();
	int pops=v[param_popSize];
	mvvv_obj.resize(maxNumrun);
	mvvv_point.resize(maxNumrun);
	for(int i=0;i<maxNumrun;i++)
	{
		mvvv_obj[i].resize(pops);
		mvvv_point[i].resize(pops);
		for(int j=0;j<pops;j++)
		{
			mvvv_obj[i][j].resize(numObj);
			mvvv_point[i][j].resize(numDim);
		}
	}
	if(glob->mp_problem->isGlobalOptKnown())
		mvv_distance.resize(maxNumrun);

	/*if (Global::g_arg[param_proName] == "FUN_FreePeak_M_OnePeak") {
		FFreePeak_M_OnePeak* mp = dynamic_cast<FFreePeak_M_OnePeak*>(glob->mp_problem.get());
		m_psr = mp->getParetoRegion();
		m_fpsr.resize(maxNumrun,vector<bool>(m_psr.size(),false));
	}*/
}

void mMultiObj::initialize(Global * glob,ParamMap &v){

    if(mMultiObj::msp_perf)        return;
   
    mMultiObj::msp_perf=unique_ptr<mMultiObj>(new mMultiObj(glob,v));
}

mMultiObj* mMultiObj::getMultiObj(){
	 return mMultiObj::msp_perf.get();
}

void mMultiObj::reInitialize(Global * glob,int pops)
{
	if(mvvv_obj[0].size()!=pops)
	{
		int maxNumrun=MAX_NUM_RUN;
		int numDim=glob->mp_problem->getNumDim();
		int numObj=glob->mp_problem->getNumObj();
		for(int i=0;i<maxNumrun;i++)
		{
			mvvv_obj[i].resize(pops);
			mvvv_point[i].resize(pops);
			for(int j=0;j<pops;j++)
			{
				mvvv_obj[i][j].resize(numObj);
				mvvv_point[i][j].resize(numDim);
			}
		}
	}
}

void mMultiObj::record(Global * glob,int index,vector<double> &obj,vector<double> &point)
{
	int ID = glob->m_runId;
	for(int i=0;i<obj.size();i++)
		mvvv_obj[ID][index][i]=obj[i];
	for(int i=0;i<point.size();i++)
		mvvv_point[ID][index][i]=point[i];
	/*if (Global::g_arg[param_proName] == "FUN_FreePeak_M_OnePeak") {
		FFreePeak_M_OnePeak* mp = dynamic_cast<FFreePeak_M_OnePeak*>(glob->mp_problem.get());
		int type = mp->getType();
		for (int i = 0; i < m_psr.size(); ++i) {
			int tidx = mp->getTreeRoot(point);
			int gidx = mp->getTree(tidx)->get_regionIdx(point);
			int bidx = mp->boxIdx(tidx, gidx);

			if (bidx == m_psr[i]) {
				if (0 == type || 1 == type)	m_fpsr[ID][i] = true;
				else  if (2 == type) {//Countable case
					MyVector v = mp->getPeak(bidx, 0).m_obj;
					if (v.getDis(obj) < 1.0e-2) {
						m_fpsr[ID][i] = true;
					}
				}
			}
		}

	}*/
}


void mMultiObj::setFileName(ParamMap &v){
	m_fileName.str("");
	
	for(auto &i:v){
		for(auto &j:Global::msm_param){
			if(i.first==param_gOptFlag||i.first==param_algId||i.first==param_proId||i.first==param_flagNoise||\
				i.first==param_flagNumPeakChange||i.first==param_flagTimeLinkage||i.first==param_numRun||\
				i.first==param_numTask||i.first==param_minNumPopSize||i.first==param_hibernatingRadius||\
				i.first==param_solutionValidationMode||i.first==param_evalCountFlag||\
				i.first==param_workingDir||i.first==param_sampleFre||i.first==param_maxEvals||i.first==param_flagNumPeakChange||\
				i.first == param_peakNumChangeMode || i.first == param_dataDirectory1) continue;
			if(i.first==j.second){			
				m_fileName<<j.first.substr(6)<<i.second<<"_";			
				break;
			}
		}		
	}
}


void mMultiObj::outputResult()
{
	stringstream ss;
	ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"PF.txt";
    ofstream out(ss.str().c_str());
	for(int i=0;i<mvvv_obj.size();i++)
	{
		for(int j=0;j<mvvv_obj[0].size();j++)
		{
			for(int z=0;z<mvvv_obj[0][0].size();z++)
				out<<mvvv_obj[i][j][z]<<" ";
			out<<endl;
		}
		out<<endl;
		out<<"*******************************************"<<endl;
	}
	out.close();
	out.clear();
	ss.str("");
	ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"PS.txt";
	out.open(ss.str().c_str());
	for(int i=0;i<mvvv_point.size();i++)
	{
		for(int j=0;j<mvvv_point[0].size();j++)
		{
			for(int z=0;z<mvvv_point[0][0].size();z++)
				out<<mvvv_point[i][j][z]<<" ";
			out<<endl;
		}
		out<<endl;
		out<<"*******************************************"<<endl;
	}
	out.close();
	out.clear();
	ss.str("");

	if(mvv_distance.size()>0)
	{
		ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"distance.txt";
		out.open(ss.str().c_str());
		vector<double> sum;
		int max=mvv_distance[0].size();
		for(int i=1;i<mvv_distance.size();i++)
			if(max<mvv_distance[i].size())
				max=mvv_distance[i].size();
		for(int i=0;i<mvv_distance.size();i++)
		{
			if(mvv_distance[i].size()<max)
			{
				int size=mvv_distance[i].size();
				double evals=mvv_distance[i][size-2];
				double temp=mvv_distance[i][size-1];
				for(int j=size;j<max;j+=2)
				{
					mvv_distance[i].push_back(evals);
					mvv_distance[i].push_back(temp);
				}
			}
		}
		sum.resize(max,0);
		for(int i=0;i<max-1;i+=2)
		{
			for(int j=0;j<mvv_distance.size();j++)
			{
				sum[i]+=mvv_distance[j][i];
				sum[i+1]+=mvv_distance[j][i+1];
			}
			sum[i]/=mvv_distance.size();
			sum[i+1]/=mvv_distance.size();
		}
		for(int i=0;i<max-1;i+=2)
			out<<sum[i]<<" "<<sum[i+1]<<endl;
		out.close();
		out.clear();
	}

	ss.str("");
	out.clear();
	ss << Global::g_arg[param_workingDir] << "Result/" << m_fileName.str() << "suc.txt";
	out.open(ss.str().c_str());
	int maxNumrun = MAX_NUM_RUN;
	double peakRate = 0;
	for (int i = 0; i < maxNumrun; ++i) {
		int j;
		for (j = 0; j < m_fpsr[i].size(); ++j) {
			if (m_fpsr[i][j] == false) break;
		}
		if (j == m_fpsr[i].size()) m_suc++;
		double rate = 0;
		for (j = 0; j < m_fpsr[i].size(); ++j) {
			if (m_fpsr[i][j] == true) rate += 1;
		}
		peakRate += rate / m_fpsr[i].size();
	}
	out << m_suc*1.0 / maxNumrun<<" "<< peakRate/ maxNumrun;
	out.close();
	
}

void mMultiObj::deleteMultiObj(){
	msp_perf.reset();
}