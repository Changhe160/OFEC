#include "LKH_OBJ.h"
#include "../LKH.h"
#include "../../../../Measure/mSingleObj.h"

unique_ptr<mLKHObj> mLKHObj::msp_perf(nullptr);

mLKHObj::mLKHObj(ParamMap &v)
{
	setFileName(v);
}

void mLKHObj::initialize(ParamMap &v)
{
	if(mLKHObj::msp_perf.get()) return;
	mLKHObj::msp_perf.reset(new mLKHObj(v));
}

void mLKHObj::setFileName(ParamMap &v){
	m_fileName.str("");
	
	for(auto &i:v){
		for(auto &j:Global::msm_param){
			if(i.first==param_gOptFlag||i.first==param_algId||i.first==param_proId||i.first==param_flagNoise||\
				i.first==param_flagNumPeakChange||i.first==param_flagTimeLinkage||i.first==param_numRun||\
				i.first==param_numTask||i.first==param_minNumPopSize||i.first==param_hibernatingRadius||\
				i.first==param_solutionValidationMode||i.first==param_evalCountFlag||\
				i.first==param_workingDir||i.first==param_sampleFre||i.first==param_maxEvals) continue;
			if(i.first==j.second){			
				m_fileName<<j.first.substr(6)<<i.second<<"_";			
				break;
			}
		}		
	}
}

mLKHObj* mLKHObj::getLKHOBJ()
{
	return mLKHObj::msp_perf.get();
}

void mLKHObj::deleteLKHOBJ()
{
	mLKHObj::msp_perf.reset();
}

void mLKHObj::outputResult()
{
	double costMin, costAve=0, costMax;
	int success=0;
	costMin=LKH::LKHAlg::mv_cost[0];
	costMax=LKH::LKHAlg::mv_cost[0];
	for(int i=0;i<LKH::LKHAlg::mv_cost.size();i++)
	{
		if(costMin>LKH::LKHAlg::mv_cost[i])
			costMin=LKH::LKHAlg::mv_cost[i];
		if(costMax<LKH::LKHAlg::mv_cost[i])
			costMax=LKH::LKHAlg::mv_cost[i];
		costAve+=LKH::LKHAlg::mv_cost[i];
		if(LKH::LKHAlg::mv_cost[i]==mSingleObj::getSingleObj()->getGOpt()[0][0])
			++success;
	}
	ostringstream os;
	os<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Sta.txt";
	ofstream out(os.str());
	out<<"Algorithm: "<<gGetAlgorithmName(Global::ms_curAlgId)<<endl;
	out<<"Number of runs: "<<LKH::LKHAlg::mv_cost.size()<<endl;
    out<<"Problem: "<<gGetProblemName(Global::ms_curProId)<<endl;
	out<<"Best: "<<costMin<<" Ave: "<<costAve/LKH::LKHAlg::mv_cost.size()<<" Worst: "<<costMax<<endl;
	out<<"Gap.Max: "<<fabs(costMax-mSingleObj::getSingleObj()->getGOpt()[0][0]);
	out<<" Gap.Ave: "<<fabs(costAve/LKH::LKHAlg::mv_cost.size()-mSingleObj::getSingleObj()->getGOpt()[0][0]);
	out<<" Gap.Min: "<<fabs(costMin-mSingleObj::getSingleObj()->getGOpt()[0][0])<<endl;
	out<<"Successes/Runs: "<<success<<"/"<<LKH::LKHAlg::mv_cost.size();
	out.close();
	out.clear();
}
