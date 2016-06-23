#ifndef OPTIMUM_H
#define OPTIMUM_H

#include "../Global/solution.h"
#include "../Utility/TypeList/TypeManip.h"

template<typename ED=CodeVReal, typename T=Solution<ED> > //T¡¡must be one derived class from Solution
class Optima{
private:
	bool m_flagLoc,m_flagGloObj; // note: m_flagGloObj only for global optimum
	vector<bool> m_foundFlag;
	vector<T> m_opt;	// note: the first optimum must be the global optimum for single obj prolems
	ProgramMode m_mode;
	bool m_recordFlag;

public:
	Optima(int dim=0,int numObj=0,bool flagLoc=false,bool flagGObj=false,ProgramMode mode=Program_Problem);
	unsigned getNumOpt()const;
	void setNumOpts(int num);
	void setLocation(int i,const vector<double> &pos);
	void setLocation(int i,const vector<TypeVar> &pos);
	void setGloObj(const vector<vector<double>> &val);
	const vector<double>& getGloObj(int i=-1)const;// i=-1 for single obj problem
	void appendOptima(const T &opt,bool flag);
	bool & flagGloObj();
	bool & flagLoc();
	void setFlagLocTrue();
	T & operator[](const unsigned id );
	const T & operator[](const unsigned id )const;
	void clear();
	void resizeDim(const int dim);
	void resizeObj(const int num);
	template<typename TS>
	bool isFound(const TS &sol,double dis,double objAcc);
	int getNumGOptFound()const;
	void printGOpt(ofstream & out)const;
	void setRecordFlag(bool f);
	void erase(int idx);
	bool isAllFound()const;
	const Solution<ED> * getNearestGOpt(const Solution<ED>& p,int *idx,double*dis)const;
	void setFoundFlagFalse(int idx=-1);
	void setFoundFlagTrue(int idx = -1);
	bool isFound(const vector<double> &obj);
};

template<typename ED, typename T>
Optima<ED,T>::Optima(int dim,int numObj,bool flagLoc,bool flagObj,ProgramMode mode):m_flagLoc(flagLoc),m_flagGloObj(flagObj),m_mode(mode),m_recordFlag(true){	
	if(m_mode==Program_Algorithm){
		m_flagLoc=false;
		m_flagGloObj=false;
		return;
	}
	if(dim>0&&numObj>0){	
		bool check[Loki::SuperSubclass<Solution<ED>,T>::value]; // check if T is derived from Solution or is type Solution
		m_opt.push_back(move(T(dim,numObj)));
		m_foundFlag.push_back(false);	
	}
}

template<typename ED, typename T>
unsigned Optima<ED,T>::getNumOpt()const{
	return m_opt.size();
}

template<typename ED, typename T>
void Optima<ED,T>::setNumOpts(int num){
	m_foundFlag.resize(num,false);
	m_opt.resize(num,m_opt[0]);
}

template<typename ED, typename T>
void Optima<ED,T>::setLocation(int i,const vector<double> &pos){
	copy(pos.begin(), pos.end(),m_opt[i].data().m_x.begin());
}


template<typename ED, typename T>
void Optima<ED,T>::setLocation(int i,const vector<TypeVar> &pos){
	m_opt[i].m_x=pos;
}

template<typename ED, typename T>
void Optima<ED,T>::setGloObj(const vector<vector<double>> &val){
	if(val.size()!=m_opt.size()) throw myException("@Optima<ED,T>::setGloObj()");
	for(int i=0;i<val.size();++i){
		copy(val[i].begin(),val[i].end(),m_opt[i].data().m_obj.begin());
	}
}

template<typename ED, typename T>
const vector<double>& Optima<ED,T>::getGloObj(int i)const{
	if(i==-1&&m_opt[0].data().m_obj.size()==1){
		if(m_opt[0].data().m_obj.size()>1){
			cout<<"warning the first one in PF was returned"<<endl;
		}
		return m_opt[0].data().m_obj;
	}else{
		return m_opt[i].data().m_obj;
	}
}

template<typename ED, typename T>
void Optima<ED,T>::appendOptima(const T &opt,bool flag){
	m_opt.push_back(move(opt));
	m_foundFlag.push_back(flag);
}

template<typename ED, typename T>
T & Optima<ED,T>::operator[](const unsigned id ){	
	return m_opt[id];
}

template<typename ED, typename T>
const T & Optima<ED,T>::operator[](const unsigned id )const{
	return m_opt[id];
}

template<typename ED, typename T>
void Optima<ED,T>::clear(){
	m_opt.clear();
	m_foundFlag.clear();
}

template<typename ED, typename T>
void Optima<ED,T>::resizeDim(const int dim){
	for(auto& i:m_opt){
		i.data().m_x.resize(dim);
	}
}
template<typename ED, typename T>
void Optima<ED, T>::resizeObj(const int num){
	for (auto& i : m_opt){
		i.data().m_obj.resize(num);
	}
}
template<typename ED, typename T>
template<typename TS>
bool Optima<ED,T>::isFound(const TS &sol,double disAcc,double objAcc){
	if(!m_recordFlag) return false;
	bool flag=true;
	if(m_mode==Program_Problem){
		if(m_flagGloObj&&!m_flagLoc){ // case: only gopt function value is provided, loc of gopt isnot provided				
			if(m_opt[0].getObjDistance_(sol.m_obj)>objAcc){
				flag=false;		
			}
		}else if(m_flagGloObj&&m_flagLoc){// case: both gopt and lopt are provided
			vector<double> dis;
			for(auto&i:m_opt) dis.push_back(i.getDistance(sol));
			if(m_opt[std::min_element(dis.begin(),dis.end())-dis.begin()].getObjDistance_(sol.m_obj)>objAcc) flag=false;
		}else{
			flag=false;
		}
	}

	if(flag){
		if(m_flagLoc){
			for(decltype (m_opt.size()) i=0;i<m_opt.size();++i){
				if(m_foundFlag[i]) continue;
				double d=m_opt[i].getDistance(sol);
				if(d<=disAcc){
					m_foundFlag[i]=true;
					return true;
				}
			}
			return false;
		}else{		
			int numFound=getNumGOptFound();
			for(int i=0;i<numFound;++i){
				if(m_opt[i].getDistance(sol)<=disAcc){
					flag=false;  // peak i has already been found
					break;
				}
			}
			if(flag){
				if(m_mode==Program_Algorithm){
					m_opt.push_back(move(T(sol)));
					m_foundFlag.push_back(true);			
				}else{
					if(numFound<m_opt.size()){
						copy(sol.m_x.begin(),sol.m_x.end(),m_opt[numFound].data().m_x.begin());
						m_foundFlag[numFound]=true;
					}else{
						cout<<"please increase distance accuracy to the global optima"<<endl;
					}
				}
			}
		}
	}
	return flag;
}

template<typename ED, typename T>
bool Optima<ED, T>::isFound(const vector<double> &obj) {
	if (!m_recordFlag) return false;	
	
	if (m_flagGloObj) { 				
		if (m_opt[0].getObjDistance_(obj)==0) {
			setFoundFlagTrue();
			return true;
		}
	}	
		
	return false;
}
template<typename ED, typename T>
inline int Optima<ED,T>::getNumGOptFound()const{
	int num=0;
	for(auto i=m_foundFlag.begin();i!=m_foundFlag.end();++i) 	if(*i==true) ++num;
	return num;
}

template<typename ED, typename T>
const Solution<ED> * Optima<ED,T>::getNearestGOpt(const Solution<ED>& p,int *idx,double*dis)const{
	
	double mindis=numeric_limits<double>::max();
	int minidx=-1;
	for(int i=0;i<m_opt.size();++i){
		if(!m_foundFlag[i]) continue;
		double dis=m_opt[i].getDistance(p);
		if(dis<mindis){
			mindis=dis;
			minidx=i;
		}
	}
	if(minidx!=-1){
		if(idx) *idx=minidx;
		if(dis) *dis=mindis;
		return &m_opt[minidx];
	}	else 	return 0;
}



template<typename ED, typename T>
void Optima<ED,T>::printGOpt(ofstream & out) const{
	for(int i=0;i<getNumGOptFound();++i){
		m_opt[i].printToFile(out);
	}
}

template<typename ED, typename T>
void Optima<ED,T>::setRecordFlag(bool f){
	m_recordFlag=f;
}
template<typename ED, typename T>
void Optima<ED,T>::erase(int idx){
	if(idx>=0&&idx<m_opt.size())	m_opt.erase(m_opt.begin()+idx); 
}
template<typename ED, typename T>
bool Optima<ED,T>::isAllFound()const{
	if(getNumGOptFound()==m_opt.size()) return true;
	return false;
}
 
template<typename ED, typename T>
bool & Optima<ED,T>::flagGloObj(){
	return m_flagGloObj;
}
template<typename ED, typename T>
bool & Optima<ED,T>::flagLoc(){
	return m_flagLoc;
}
template<typename ED, typename T>
void Optima<ED,T>::setFlagLocTrue(){
	m_flagGloObj=true;
	m_flagLoc=true;
}

template<typename ED, typename T>
void Optima<ED, T>::setFoundFlagFalse(int idx){
	if (idx == -1) for (size_t i = 0; i != m_foundFlag.size(); ++i)  m_foundFlag[i] = false;
	else m_foundFlag[idx] = false;
}

template<typename ED, typename T>
void Optima<ED, T>::setFoundFlagTrue(int idx){
	if (idx == -1) for (size_t i = 0; i != m_foundFlag.size(); ++i)  m_foundFlag[i] = true;
	else m_foundFlag[idx] = true;
}
#endif