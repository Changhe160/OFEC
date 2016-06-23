#ifndef ENCONDING_H
#define ENCONDING_H

#include "../Utility/include.h"
#include "../Utility/TypeVar/typeVar.h"
#include "../Utility/TypeList/Typelist.h"

enum Encoding{Code_VReal=0,Code_VInt,Code_VBinary,Code_VVar};

struct VirtualEncoding{
	vector<double> m_obj;
	VirtualEncoding(size_t s=0):m_obj(s){}
	VirtualEncoding(const vector<double> &o) :m_obj(o){}
	virtual ~VirtualEncoding()=0;
	double getObjDistance(const vector<double>&rhs, DistanceMode mode = DIS_MANHATTAN)const{
		if (this->m_obj.size() != rhs.size()){
			throw myException("the number of objetives must be the same @ double VirtualEncoding::getObjDistance()");
		}
		double dis = 0;
		switch (mode){
			case DIS_EUCLIDEAN:
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					dis += (this->m_obj[i] - rhs[i])*(this->m_obj[i] - rhs[i]);
				}
				dis=sqrt(dis);
				break;
			case DIS_MANHATTAN:
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					dis += fabs(this->m_obj[i] - rhs[i]);
				}
				break;			
			case DIS_HAMMING:			
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					if(this->m_obj[i] !=rhs[i]) dis+=1;
				}
				break;	
		}		
		return dis;
	}
};

inline VirtualEncoding::~VirtualEncoding(){ }

struct CodeVReal:public VirtualEncoding
{
	vector<double> m_x;
	CodeVReal() = default;
	CodeVReal(size_t s1,size_t s2):VirtualEncoding(s2),m_x(s1){}
	CodeVReal(int s1,int s2):VirtualEncoding(s2),m_x(s1){}
	CodeVReal(const vector<double> &x, const vector<double> &o) :VirtualEncoding(o), m_x(x){}
	CodeVReal(const vector<double> &x, const int s) :VirtualEncoding(s), m_x(x){}
	template<typename Itor>
	CodeVReal(Itor begin,Itor end):m_x(begin,end){}
	double operator[]( int i)const{
		return m_x[i];
	}
	double& operator[]( int i){
		return m_x[i];
	}
};

struct CodeVInt:public VirtualEncoding
{
	vector<int> m_x;
	CodeVInt(size_t s1,size_t s2):VirtualEncoding(s2),m_x(s1){}
	CodeVInt(int s1,int s2):VirtualEncoding(s2),m_x(s1){}
	template<typename Itor >
	CodeVInt(Itor  begin,Itor end):m_x(begin,end){}
	int operator[]( int i)const{
		return m_x[i];
	}
	int& operator[]( int i){
		return m_x[i];
	}
};

struct CodeVVar:public VirtualEncoding
{
	vector<TypeVar> m_x;
	CodeVVar(size_t s1,size_t s2):VirtualEncoding(s2),m_x(s1){}
	CodeVVar(int s1,int s2):VirtualEncoding(s2),m_x(s1){}
	template<typename Itor>
	CodeVVar(Itor begin,Itor end):m_x(begin,end){}
	const TypeVar& operator[]( int i)const{
		return m_x[i];
	}
	TypeVar& operator[]( int i){
		return m_x[i];
	}
};

struct CodeVBinary:public VirtualEncoding
{
	vector<bool> m_x;
	CodeVBinary(size_t s1,size_t s2):VirtualEncoding(s2),m_x(s1){}
	CodeVBinary(int s1,int s2):VirtualEncoding(s2),m_x(s1){}
	bool operator[]( int i)const{
		return m_x[i];
	}
	struct Proxy {
		Proxy &operator=(Proxy t){
			
		}
		Proxy& operator=(bool t){
			m_ref[m_i]=t;
			return *this;
		}
		operator bool() const{
			return m_ref[m_i];
		}
		Proxy(int i,vector<bool>&x):m_i(i),m_ref(x){}
		vector<bool> &m_ref;
		int m_i;
	 };
	friend struct Proxy;
	Proxy operator[]( int i){
		return Proxy(i,m_x);
	}
};


//please add your encoding scheme below and choose it as the default encoding 
typedef LOKI_TYPELIST_4(CodeVReal,CodeVInt,CodeVBinary,CodeVVar) EncodingTypeList;

#define ENCODING_CODE Code_VReal

typedef Loki::TL::TypeAt<EncodingTypeList,ENCODING_CODE>::Result EncodingType;


#endif