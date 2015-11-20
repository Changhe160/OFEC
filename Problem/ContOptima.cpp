#include "ContOptima.h"
#include "ContinuousProblem.h"
#include "../Utility/definition.h"
extern boost::mutex g_mutex;

#ifdef OFEC_DEMON
#include "../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

ContOptimum::ContOptimum(int dim,int numObj):Solution<CodeVReal>(dim,numObj),m_radius(0),m_startRadius(0),m_vbest(dim),m_step(0),m_order(1),m_deratingObj(numObj){	

}
ContOptimum::ContOptimum(const CodeVReal & s,double startRadi,double step):Solution<CodeVReal>(s),m_radius(0),m_startRadius(startRadi),m_vbest(s.m_x),m_step(step),\
	m_isReady(false),m_order(1),m_deratingObj(s.m_obj.size()){ 
	m_deratingObj=Solution<CodeVReal>::getWorstSolutionSoFar().obj();
	if(m_startRadius==0){
		m_startRadius=CAST_PROBLEM_CONT->getDisAccuracy()/3;
	}
	if(m_step==0){
		m_step=m_startRadius/3;
	}
	generateSample();	
}
ContOptimum::ContOptimum(const Solution<CodeVReal> & s,double startRadi,double step):Solution<CodeVReal>(s),m_radius(0),m_startRadius(startRadi),m_vbest(s.data().m_x),m_step(step),\
	m_isReady(false),m_order(1),m_deratingObj(s.obj().size()){ 
	m_deratingObj=Solution<CodeVReal>::getWorstSolutionSoFar().obj();
	if(m_startRadius==0){
		m_startRadius=CAST_PROBLEM_CONT->getDisAccuracy()/3;
	}
	if(m_step==0){
		m_step=m_startRadius/3;
	}
	generateSample();	
}
void ContOptimum::setRadius(double r){
	m_radius=r;
}
double ContOptimum::getRadius(){
	return m_radius;
}

void ContOptimum::setDeratingObj(int num){
	m_order=num;
}
double ContOptimum::getRefRadi(const Solution<CodeVReal> &s){
	if(!m_isReady) return 0;
	{
		#ifdef OFEC_DEMON
		Ulock lock(m_mutex);
		#endif
		MyVector v(s.data().m_x);
		v-=m_vbest;
		if(v.length()==0) return m_radius;
		vector<double> angle;
		for(auto&i:m_sample){
			angle.push_back(v.getAngle(i));
		}
		int idx=min_element(angle.begin(),angle.end())-angle.begin();
		double mang=angle.empty()?OFEC_PI:angle[idx];
		if(mang>3*OFEC_PI/180){		
			m_isReady=false;
			v.normalize();
			{
				Ulock lock2(g_mutex);
				creatOneSample(v);
			}
			m_isReady=true;
			return m_sample.back().length();	
		}
		return m_sample[idx].length();
	}
}

void ContOptimum::generateSample(){
	{
		#ifdef OFEC_DEMON
		Ulock lock(m_mutex);
		Ulock lock2(g_mutex);
		#endif
		for(int tr=0;tr<GET_NUM_DIM*2;++tr){
			MyVector vnor(GET_NUM_DIM);
			vnor.randomize(-1,1);
			vnor.normalize();		
			creatOneSample(vnor);		
		}	
		m_isReady=true;
	}
	#ifdef OFEC_DEMON
	msp_buffer->updateFitnessLandsacpe_();
	#endif
}

void ContOptimum::creatOneSample(const MyVector &vnor){
	
	MyVector vtr(vnor);
	Solution s0(*this),s1(*this);
	double r=m_startRadius;
	bool flag=true;
	do{
		s0=s1;
		vtr=vnor*r;
		vtr+=m_vbest;
		copy(vtr.begin(),vtr.end(),s1.data().m_x.begin());
		if(Global::msp_global->mp_problem->isValid(s1)){
			s1.evaluate(false);
			Global::msp_global->mp_problem->cevals()++;
		}else{	 	
			SolutionValidation mode=VALIDATION_SETTOBOUND;
			s1.validate(&mode);
			s0=s1;
			flag=false;
			break;	
		}			
		r+=m_step;
	}while(s0>=s1);
	if(flag){
		binarySearch(s0,s1);
	}
	MyVector vr(s0.data().m_x);
	vr-=m_vbest;
	m_sample.push_back(move(vr));
	if(m_sample.size()==1||m_radius>r) m_radius=r;
}

ContOptimum& ContOptimum::operator=(const ContOptimum& rhs){
	Solution::operator=(rhs);
	m_radius=rhs.m_radius;
	m_sample=rhs.m_sample;
	m_startRadius=rhs.m_startRadius;
	m_vbest=rhs.m_vbest;
	m_deratingObj=rhs.m_deratingObj;
	m_step=rhs.m_step;
	m_isReady=rhs.m_isReady;
	m_order=rhs.m_order;
	return *this;
}

ContOptimum::ContOptimum(const ContOptimum& rhs):Solution(rhs){
	m_radius=rhs.m_radius;
	m_sample=rhs.m_sample;
	m_startRadius=rhs.m_startRadius;
	m_vbest=rhs.m_vbest;
	m_deratingObj=rhs.m_deratingObj;
	m_step=rhs.m_step;
	m_isReady=rhs.m_isReady;
	m_order=rhs.m_order;
}
MyVector& ContOptimum::operator[](int i){
	return m_sample[i];
}
int ContOptimum::size(){
	return m_sample.size();
}
MyVector& ContOptimum::getBest(){
	return m_vbest;
}
void ContOptimum::binarySearch(Solution<CodeVReal> &s0, Solution<CodeVReal> &s1){
	Solution<CodeVReal> s;
	for(int i=0;i<GET_NUM_DIM;++i){
		s.data().m_x[i]=(s0.data().m_x[i]+s1.data().m_x[i])/2.;
	}
	s.evaluate(false);
	if(s0.getDistance(s)<=1.e-3){
		s0=s;
		return;
	}
	if(s<s0) binarySearch(s,s1);
	if(s>s0) binarySearch(s0,s);
}
