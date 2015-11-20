#ifndef FAMFRAMEWORK_POP_H
#define FAMFRAMEWORK_POP_H
#include "FAMFDera.h"

template<typename, typename, typename> class FAMF;
class FAMFPopPSO;
template<typename TypeIndi,typename TypeMainPop>
class FAMFPop: public TypeMainPop{	
	template<typename, typename, typename> friend class FAMF;
	friend class FAMFPopPSO;
	enum{f_hiber=0};
protected:
	int m_countRadiusDecreace;
	bool m_stagnantFlag;
public:	
	template<typename ED>
	FAMFPop(Group<ED,TypeIndi> &g);
	bool isStagnant(double degree,double avgRadius);
	void estimateInitialRadius();
	virtual ~FAMFPop(){}
	void updateCurRadius(bool mode=false);
	bool updateHiberFlag( );
	bool isHibernating();
	void wakeup();
	void reduceRadius(double agvr);
	virtual void deratingFitness(unique_ptr<TypeMainPop>& single);
	void computeCenter();
	bool isStagnant();
	void setHiberFlag(bool flag);
};

template<typename TypeIndi,typename TypeMainPop>
template<typename ED>
FAMFPop<TypeIndi,TypeMainPop>::FAMFPop( Group<ED,TypeIndi> &gr):TypeMainPop(gr),m_countRadiusDecreace(0),m_stagnantFlag(false){
	estimateInitialRadius();
}

template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::estimateInitialRadius(){
	if(this->m_popsize<=1) {
		this->m_initialRadius=0;
		return;
	}
	this->updateCurRadius(true);
	this->m_initialRadius=this->m_curRadius;
}

template<typename TypeIndi,typename TypeMainPop>
bool FAMFPop<TypeIndi,TypeMainPop>::isStagnant(double degree,double avgRadius){
	if(m_countRadiusDecreace>=1.*this->m_popsize&&this->m_curRadius>=avgRadius&&this->m_curRadius>degree*CAST_PROBLEM_CONT->getDomainSize()||m_countRadiusDecreace>=10*this->m_popsize) m_stagnantFlag= true;
	else m_stagnantFlag= false;
	return m_stagnantFlag;
}

template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::updateCurRadius(bool mode){
		//mode=true for distance between representatice to the center; else between pself to the center		
		TypeMainPop::updateCurRadius(mode);
		if(this->m_popsize<2) return;
		if(Global::ms_curProId==Global::msm_pro["DYN_CONT_MovingPeak"]){			
			double qua=this->m_radiusQaulity;
			this->evaluateRadius(this->m_curRadius);
			if(qua>this->m_radiusQaulity) this->m_radiusQaulity=qua;
			else this->m_initialRadius=this->m_curRadius;			
		}
}

template<typename TypeIndi,typename TypeMainPop>
bool FAMFPop<TypeIndi,TypeMainPop>::updateHiberFlag( ){
	if(this->m_flag[f_hiber]) return false;
	int idx=-1;
	if(CAST_PROBLEM_CONT->getGOpt().flagLoc()){	
		CAST_PROBLEM_CONT->getGOpt().getNearestGOpt(*this->m_best[0],&idx,0);
	}else	if(CAST_PROBLEM_CONT->getGOpt().flagGloObj()){
		idx=0;
	}	
	if(idx!=-1){		
		double dis=this->m_best[0]->getObjDistance_(CAST_PROBLEM_CONT->getGOpt().getGloObj(idx));
		double disx=0;
		if(CAST_PROBLEM_CONT->getGOpt().flagLoc()) disx=this->m_best[0]->getDistance(CAST_PROBLEM_CONT->getGOpt()[idx]);
		if(disx<=CAST_PROBLEM_CONT->getDisAccuracy()&&dis<=CAST_PROBLEM_CONT->getAccuracy()){
			this->m_flag[f_hiber]=true;
		}else this->m_flag[f_hiber]=false;
	}
	return this->m_flag[f_hiber];
}

template<typename TypeIndi,typename TypeMainPop>
bool FAMFPop<TypeIndi,TypeMainPop>::isHibernating(){
	return this->m_flag[f_hiber];
}

template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::wakeup(){
	this->m_flag[f_hiber]=false;
}

template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::reduceRadius(double avgr){
		this->m_initialRadius=avgr;
}
template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::deratingFitness(unique_ptr<TypeMainPop>& single){

	for(auto &i:this->m_best){
		FAMFDerating::derateFitness(*i);
	}
	FAMFDerating::derateFitness(this->m_center);

	for(auto &i:this->m_pop) {
		bool flag=FAMFDerating::derateFitness(i->self());
		if(&i->self()!=&i->representative()){
			FAMFDerating::derateFitness(i->representative());
		}
		if(!flag) {
			for(auto &best:this->m_best){
				if(i->self()>best->self()){
					single->add(i.get(),false);
					break;
				}
			}
		}
	}
	
}
template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::computeCenter(){
	if(this->m_popsize<1) return;
		
	if(Global::msp_global->mp_problem->isProTag(CONT)){
		for( int i=0;i<GET_NUM_DIM;i++){
			double x=0.;
			for(int j=0;j<this->m_popsize;j++){
				Solution<CodeVReal> & chr=this->m_pop[j]->representative();
				x=chr.data()[i]+x;
			}
			this->m_center.data()[i]=x/this->m_popsize;
		}			
		this->m_center.evaluate(false);
		
		
		FAMFDerating::derateFitness(this->m_center);
		for(auto &i:this->m_best){
			FAMFDerating::derateFitness(*i);
		}
		this->updateBestArchive(this->m_center);
	}

}

template<typename TypeIndi,typename TypeMainPop>
bool FAMFPop<TypeIndi,TypeMainPop>::isStagnant(){
	return m_stagnantFlag;
}
template<typename TypeIndi,typename TypeMainPop>
void FAMFPop<TypeIndi,TypeMainPop>::setHiberFlag(bool flag){
	this->m_flag[f_hiber]=flag;
}
#endif