#include "mROOT.h"

#include "../Problem/DOP/DynamicProblem.h"

unique_ptr<mROOT> mROOT::msp_root(nullptr);
mROOT::mROOT(ParamMap &v):m_timeWindow((int)(v[param_timeWindow])),m_numChanges((int)MAX_NUM_RUN,0),\
	mq_solution((int)MAX_NUM_RUN),m_mean(0),m_std(0),mv_root((int)MAX_NUM_RUN){

}

mROOT::~mROOT(){
	mv_root.clear();
	mv_avgRoot.clear();
	mq_solution.clear();
}
bool mROOT::initilizeROOT(ParamMap &v){
	msp_root.reset(new mROOT(v));
	if(msp_root.get()){	return true;
	}else return false;
}
void mROOT::deleteRoot(){
	msp_root.reset();
}
bool mROOT::updateObjValue(Global* glob){
    if(!msp_root) return false;
	vector<double> obj;
	glob->mp_problem->getObjGlobalOpt(obj);
	int runId=glob->m_runId;
	for(unsigned i=0; i<mq_solution[runId].size();i++){
		for(unsigned j=0;j<mq_solution[runId][i].size();j++){
			if(mq_solution[runId][i][j].getNumDim()!=glob->mp_problem->getNumDim()) return false;
			vector<double> obj_=mq_solution[runId][i][j].obj();
			mq_solution[runId][i][j].evaluate(false);
			obj_.push_back(fabs(obj[0]-mq_solution[runId][i][j].obj(0)));
			mq_solution[runId][i][j].obj()=obj_;
		}
	}
	return true;
}
double mROOT::getMeanRoot(){

	return m_mean;
}
void mROOT::calculateMean(){

	for(unsigned runId=0;runId<mv_root.size();runId++){
	//compute the root value of previous run
		double sum=0;
		for(unsigned i=0;i<mv_root[runId].size();i++){
			sum+=mv_root[runId][i];
		}
		mv_avgRoot.push_back(sum/mv_root[runId].size());
	}

	m_mean=0;
	for(unsigned i=0;i<mv_avgRoot.size();i++){
		m_mean+=mv_avgRoot[i];
	}
	m_mean/=mv_avgRoot.size();
	m_std=0;
	for(unsigned i=0;i<mv_avgRoot.size();i++){
		m_std+=(mv_avgRoot[i]-m_mean)*(mv_avgRoot[i]-m_mean);
	}
	m_std=sqrt(m_std/mv_avgRoot.size());
}
double mROOT::getROOTStd(){

	return m_std;
}
void mROOT::record(const Solution<> & chr,Global *glob){
    if(!msp_root) return;
	int runId=glob->m_runId;

	int changes=CAST_PROBLEM_DYN->getChangeCounter();
	if(mq_solution[runId].size()>=m_timeWindow&&m_numChanges[runId]<=changes){
		// calculate the root value of solutions going to be poped out
		double minSum=numeric_limits<double>::max();
		for(unsigned i=0;i<mq_solution[runId].front().size();i++){
			double sum=0;
			for(unsigned j=0;j<mq_solution[runId].front()[i].obj().size();j++){
				sum+=mq_solution[runId].front()[i].obj()[j];
			}
			if(minSum>sum) minSum=sum;
		}
		mv_root[runId].push_back(minSum/m_timeWindow);
		mq_solution[runId].pop_front();
	}

	if(m_numChanges[runId]<=changes){
		// push a new element into the queue
		mq_solution[runId].push_back(move(vector<Solution<>>()));
		m_numChanges[runId]++;
	}
	mq_solution[runId].back().push_back(chr);
	vector<double> obj;
	glob->mp_problem->getObjGlobalOpt(obj);
	vector<double> obj_=mq_solution[runId][mq_solution[runId].size()-1][mq_solution[runId][mq_solution[runId].size()-1].size()-1].obj();
	mq_solution[runId][mq_solution[runId].size()-1][mq_solution[runId][mq_solution[runId].size()-1].size()-1].obj()[0]=fabs(obj_[0]-obj[0]);
	
}
