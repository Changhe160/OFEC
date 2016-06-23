#include "test.h"
#include "initialization.h"
#include "run.h"
#include "../Problem/Combination/TSP/OptimalEdgeInfo.h"
#include "../Algorithm/Other/LKH/LKH_outPutResult/LKH_OBJ.h"

mutex g_mutex;
mutex g_mutexStream;
void outputResults(){
#ifndef OFEC_PROBLEM_DEBUG
	if(mLKHObj::getLKHOBJ()!=nullptr)
	{
		mLKHObj::getLKHOBJ()->outputResult();
		mLKHObj::deleteLKHOBJ();
		return ;
	}

    if(mSingleObj::getSingleObj()!=nullptr){
            mSingleObj::getSingleObj()->outputResult();
			if(OptimalEdgeInfo::getOptimalEdgeInfo()!=nullptr){
				OptimalEdgeInfo::getOptimalEdgeInfo()->output();
				OptimalEdgeInfo::deleteOptimalEdgeInfo();
			}
			mSingleObj::deleteSingleObj();
    }

    if(mMultiModal::getPopInfor()!=nullptr){
            mMultiModal::getPopInfor()->output();
			mMultiModal::deleteMultiModal();
    }

	if(mMultiObj::getMultiObj()!=nullptr){
            mMultiObj::getMultiObj()->outputResult();
			mMultiObj::deleteMultiObj();
    }
#endif
}

int go(vector<int> runId){

	for(auto & r:runId){
		try {
			
			g_mutex.lock();
			Run run(r,Global::g_arg);
			g_mutex.unlock();
		
		//run.test();
			run.go();
			/*g_mutexStream.lock();
			cout << "run " << r << " finishes "<< endl;
			g_mutexStream.unlock();*/
		}
		catch(exception &e){
			g_mutexStream.lock();
			cout<<"exception, runId "<<r<<": "<<e.what()<<endl;
			g_mutexStream.unlock();
		}
	}
	return 0;
}

void run(){
	int numTask=0;
	if(Global::g_arg.find(param_numTask)!=Global::g_arg.end()) numTask=(Global::g_arg[param_numTask]);
	if(numTask==1)	{
		//***** non-concurrent run******//
		//cout<<"Warning: Program is running in SINGLE-THREAD mode!"<<endl;
		for(auto i=0;MAX_NUM_RUN>i;i++){
			vector<int> runs;
			runs.push_back(i);
			go(runs);
		}
		return;
	}else if(numTask==0){ // auto mode, the number of threads depends on the number of logical/physical threads 
		numTask = thread::hardware_concurrency();
	}

	//**** concurrent run *****//
	//cout<<"Warning: Program is running in MULTI-THREAD mode!"<<endl;
	if (MAX_NUM_RUN<numTask) numTask = MAX_NUM_RUN;

	vector<thread> atrd;
	int rest = MAX_NUM_RUN%numTask;
	int id1 = 0, id2 = id1 + MAX_NUM_RUN / numTask - 1 + (rest-->0 ? 1 : 0);
	for (int i = 0; i<numTask; i++) {
		vector<int> runs;
		for (int r = id1; r <= id2; r++) runs.push_back(r);
		id1 = id2 + 1;
		id2 = id1 + MAX_NUM_RUN / numTask - 1 + (rest-->0 ? 1 : 0);
		atrd.push_back(thread(go, runs));
	}

	for (auto&t : atrd) t.join();
	
}

int main(int argc,char* argv[]){

	time_t timer_start,timer_end;
	setGlobalParameters(argc,argv);
	registerProNdAlg();
	setAddtParameters();
	time(&timer_start);
	run();
	time(&timer_end);
	cout<<"Time used: "<<difftime(timer_end,timer_start)<<" seconds"<<endl;	
    outputResults();
	
	return 0;
}