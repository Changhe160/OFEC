#ifndef CONTINUOUSE_OPTIMA
#define CONTINUOUSE_OPTIMA
#include "../Global/solution.h"
#include "../Utility/myVector.h"
class ContOptimum: public Solution<CodeVReal>{
protected:	
	double m_radius;
	vector<MyVector> m_sample;
	double m_startRadius;
	MyVector m_vbest;
	vector<double> m_deratingObj;
	double m_step;
	bool m_isReady;
	int m_order;
#ifdef OFEC_DEMON
	boost::mutex m_mutex;
#endif // OFEC_DEMON

public:
	ContOptimum(int dim,int numObj=1);
	ContOptimum(const Solution<CodeVReal> & s,double r=0,double step=0);
	ContOptimum(const CodeVReal & s,double r=0,double step=0);
	void setRadius(double r);
	double getRefRadi(const Solution<CodeVReal> &s);
	double getRadius();
	void setDeratingObj(int num);
	vector<double> & getDeratingObj();
	ContOptimum& operator=(const ContOptimum& rhs);
	ContOptimum(const ContOptimum& rhs);
	MyVector& operator[](int i);
	int size();
	MyVector& getBest();
protected:
	void generateSample();
	void creatOneSample(const MyVector &vnor);
	void binarySearch(Solution<CodeVReal> &s0, Solution<CodeVReal> &s1);
	
};

inline vector<double> &ContOptimum::getDeratingObj(){
	double w=pow(1.001,-sqrt(m_order));
	for(int i=0;i<m_obj.size();++i){
		m_deratingObj[i]=Solution<CodeVReal>::getBestSolutionSoFar().obj(i)*(1-w)+Solution<CodeVReal>::getWorstSolutionSoFar().obj(i)*w;	
	}
	return m_deratingObj;
}
#endif