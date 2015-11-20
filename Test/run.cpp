#include "run.h"
#include "test.h"
#include "../Problem/Combination/TSP/OptimalEdgeInfo.h"
#include "../Algorithm/Other/LKH/LKH_outPutResult/LKH_OBJ.h"

Run::Run(int runId,ParamMap &v):m_runId(runId){
	if(!Global::msp_global.get()){
		#ifdef OFEC_CONSOLE
		Global::msp_global.reset(new Global(m_runId,1./7,(m_runId+1.)/((int)(MAX_NUM_RUN)+1.)));
		#endif
		#ifdef OFEC_DEMON
		Global::msp_global.reset(new Global(m_runId,1./7,1./11));
		#endif	
		Global::msp_global->mp_problem.reset(Global::ms_classFactory.constructProblem((v[param_proName]))(v));
		
		#ifdef OFEC_CONSOLE
		if(Global::msp_global->mp_problem->isProTag(SOP)){
			if(Global::msp_global->mp_problem->isProTag(CONT)){
				if(CAST_PROBLEM_CONT->getGOpt().flagGloObj()) Global::g_arg[param_gOptFlag]=true;
				else  Global::g_arg[param_gOptFlag]=false;
			}else if(Global::msp_global->mp_problem->isProTag(TSP)){
				if(CAST_TSP->getGOpt().flagGloObj()) Global::g_arg[param_gOptFlag]=true;
				else  Global::g_arg[param_gOptFlag]=false;
			}

			if(mSingleObj::getSingleObj()==nullptr){
				if(Global::msp_global->mp_problem->isProTag(DOP))
					mSingleObjDyn::initialize(Global::g_arg); 
				else
					mSingleObj::initialize(Global::g_arg); 
				mSingleObj::getSingleObj()->setFileName(Global::g_arg);
			}
			
			vector<double> gOpt;
			if(mSingleObj::getSingleObj() &&Global::msp_global->mp_problem->getObjGlobalOpt(gOpt)){
				mSingleObj::getSingleObj()->addGOpt(m_runId,gOpt[0]);
			}
			
			if(m_runId==0){
					mSingleObj::getSingleObj()->setAccuracy(Global::msp_global->mp_problem->getAccuracy());
					mSingleObj::getSingleObj()->setCompareType(Global::msp_global->mp_problem->getOptType());
					mSingleObj::getSingleObj()->setProParameter(Global::msp_global->mp_problem->m_proPar);
					
			}
			
		}
		if(Global::msp_global->mp_problem->isProTag(MMP)){
			if(mMultiModal::getPopInfor()==nullptr){
					mMultiModal::initialize(Global::g_arg);
					mMultiModal::getPopInfor()->setFileName(Global::g_arg);	
			}
		}
		if(Global::msp_global->mp_problem->isProTag(MOP)){
			if(mMultiObj::getMultiObj()==nullptr){	
				if(Global::msp_global->mp_problem->isProTag(DOP))
					mMultiObj::initialize(Global::msp_global.get(),Global::g_arg);
				else mMultiObj::initialize(Global::msp_global.get(),Global::g_arg); 
				mMultiObj::getMultiObj()->setFileName(Global::g_arg);
			}
		}
		if(Global::msp_global->mp_problem->isProTag(TSP))
		{
			if(OptimalEdgeInfo::getOptimalEdgeInfo()==nullptr)
				OptimalEdgeInfo::initialize(v);
			if(IS_ALG_NAME(Global::ms_curAlgId,"ALG_LKH"))
			{
				if(mLKHObj::getLKHOBJ()==nullptr)
					mLKHObj::initialize(v);
			}
		}
		
		#endif
		#ifndef OFEC_PROBLEM_DEBUG
			Global::msp_global->mp_algorithm.reset(Global::ms_classFactory.constructAlgorithm((v[param_algName]))(v));
		#endif
		#ifdef OFEC_CONSOLE
			if (mSingleObj::getSingleObj() && Global::msp_global->mp_algorithm)	mSingleObj::getSingleObj()->setAlgParameter(Global::msp_global->mp_algorithm->m_algPar);
		#endif
	}

}
Run::~Run(){
	Global::msp_global.reset();	
}
ReturnFlag Run::go(){
	return Global::msp_global->mp_algorithm->run();
}

void Run::test(){
	
	//FreePeak_D_M_OnePeak
	/*Solution<CodeVReal> x(Global::msp_global->mp_problem->getNumDim(), Global::msp_global->mp_problem->getNumObj());
	Global::msp_global->mp_problem->initializeSolution(x.data());
	for (int i = 0; i < 10000; ++i){
		Global::msp_global->mp_problem->evaluate(x,true,Program_Algorithm);
	}*/
	
	//NT=1 NR=1 NO=1 ND=5 NB=100 PN=FUN_FreePeak_D_OnePeak AN=ALG_FAMF_PSO PS=1 NC=1 CF=1 SL=1.0 FNPC=false PPB=1 NGO=1 TW=1 CT=1 CT2=0 PF=false NF=false CR=1.0 CoF=0.005
	//NT=1 NR=1 ND=5 NP=100 PN=DYN_CONT_MovingPeak AN=ALG_FAMF_PSO NC=1 SL=1.0 CT=0 CR=1.0 PS=1 CoF=0.005 CF=5000 PIM=0 RBP=3 
	//NT=1 NR=1 ND=5 NP=100 PN=DYN_CONT_RotationDBG AN=ALG_FAMF_PSO NC=1 CT=0 CR=1.0 PS=1 CoF=0.005 CF=5000 PIM=0 RBP=3 
	//FreePeak_D_OnePeak, MovingPeak,RotationDBG,ComPositionDBG
	Solution<CodeVReal> x(Global::msp_global->mp_problem->getNumDim(), Global::msp_global->mp_problem->getNumObj());
	auto timer_start=std::chrono::system_clock::now();
	for (int i = 0; i < 1000000; ++i){
		Global::msp_global->mp_problem->initializeSolution(x.data());
		Global::msp_global->mp_problem->evaluate(x, false, Program_Problem);
	}
	auto timer_end = std::chrono::system_clock::now();
	auto td = timer_end - timer_start;
	
	string ss = Global::g_arg[param_workingDir];
	ss += "Result/";
	ss += mSingleObj::getSingleObj()->m_fileName.str();
	ss += "time.txt";
	static long long atd = 0;
	static float invisble = 0;
	float num = 0;
	if (Global::g_arg[param_numDim] == 5){
		if (Global::g_arg[param_proName] == "DYN_CONT_MovingPeak" || Global::g_arg[param_proName] == "DYN_CONT_RotationDBG"){
			for (int i = 0; i < 1000; ++i){
				DynamicContinuous*pro = dynamic_cast<DynamicContinuous*>(Global::msp_global->mp_problem.get());
				pro->change();
				num += pro->getNumberofPeak() - pro->getNumofVisablePeaks();
			}
			num /= 1000;
		}
		else if (Global::g_arg[param_proName] == "FreePeak_D_OnePeak"){

		}
	}

	g_mutex.lock();
	ofstream out(ss.c_str());
	atd+=std::chrono::duration_cast<std::chrono::milliseconds>(td).count();
	invisble += num;
	//cout << "time used (milliseciond): " << std::chrono::duration_cast<std::chrono::milliseconds>(td).count() << endl;
	out << atd/(int)Global::g_arg[param_numRun] << endl;
	out << invisble / (int)Global::g_arg[param_numRun] << endl;
	out.close();
	g_mutex.unlock();
	
}