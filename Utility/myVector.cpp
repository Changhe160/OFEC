#include "myVector.h"
#include "myexcept.h"
#include "../Global/global.h"

MyVector::MyVectorProxy::MyVectorProxy(MyVector& v,int idx):m_vec(v),m_idx(idx){

}
MyVector::MyVectorProxy& MyVector::MyVectorProxy::operator=(const MyVectorProxy& rhs){
	m_vec.m_data[m_idx]=rhs.m_vec.m_data[rhs.m_idx];
	m_vec.m_write=true;
	return *this;
}
MyVector::MyVectorProxy::operator double ()const{
	return m_vec.m_data[m_idx];
}
MyVector::MyVectorProxy&  MyVector::MyVectorProxy::operator=(double val){
	m_vec.m_data[m_idx]=val;
	m_vec.m_write=true;
	return *this;
}

MyVector::MyVectorProxy & MyVector::MyVectorProxy::operator+=(double val){
	m_vec.m_data[m_idx]+=val;
	m_vec.m_write=true;
	return *this;
}

MyVector::MyVector(const vector<TypeVar> &v):m_length(0),m_write(true){
	for(auto&i:v) {
		if(i.is_double()){ 
			m_data.push_back(double(i));		
		}else{
			throw myException("cannot assign non-double tpyeVar to double@MyVector(const Solution &s)");
		}
	}
}
MyVector::MyVector( int size,double *val):m_data(size),m_length(0),m_write(true){
	copy(val,val+size,m_data.begin());
}
MyVector::MyVector( int size):m_data(size),m_length(0),m_write(true){

}
MyVector::MyVector( int size,double val):m_data(size,val),m_length(0),m_write(true){
	
}
MyVector::MyVector(const vector<double> & v):m_data(v),m_length(0),m_write(true){

}
MyVector::~MyVector(){
	m_data.clear();
}
MyVector::MyVectorProxy MyVector::operator [](int idx) {
	return MyVector::MyVectorProxy(*this,idx);
}
const double & MyVector::operator [](int idx)const{
	return m_data[idx];
}
MyVector &MyVector::operator =(const MyVector & v){
	if(this==&v) return *this;
	m_data=v.m_data;
	m_length=v.m_length;
	m_write=v.m_write;
	return *this;
}
MyVector &  MyVector::operator +=(const MyVector & v){
	if(m_data.size()!=v.m_data.size()) throw myException("the size of two vectors must be same by + operation@MyVector::operator +=");
	transform(m_data.begin(),m_data.end(),v.m_data.begin(),m_data.begin(),plus<double>());
	m_write=true;
	return *this;
}
 MyVector &  MyVector::operator -=(const MyVector & v){
	if(m_data.size()!=v.m_data.size()) throw myException("the size of two vectors must be same by - operation@MyVector::operator -=");
	transform(m_data.begin(),m_data.end(),v.m_data.begin(),m_data.begin(),minus<double>());
	m_write=true;
	return *this;
}
 double MyVector::operator *(const MyVector & v){
	if(m_data.size()!=v.m_data.size()) throw myException("the size of two vectors must be same by * operation@MyVector::operator *");
	double sum=0;
	for(unsigned i=0;i<m_data.size();i++) sum+=m_data[i]*v.m_data[i];
	return sum;
}
 void MyVector::projectionToV(const MyVector &v){

	double sv=0,vv=0;
	for(unsigned i=0;i<m_data.size();i++){
		sv+=m_data[i]*v.m_data[i];
		vv+=v.m_data[i]*v.m_data[i];
	}
	for(unsigned i=0;i<m_data.size();i++){
		m_data[i]=sv*v.m_data[i]/vv;
	}
	m_write=true;
}

 void MyVector::normalize(){
	if (length() == 0) throw myException("vector length is zero @ MyVector::normalize()");
	double sum=0;
	for_each(m_data.begin(),m_data.end(),[&sum](double v){sum+=v*v;});
	sum=sqrt(sum);
	for_each(m_data.begin(),m_data.end(),[&sum](double &v){v=v/sum;});
	m_length=1;
	m_write=false;
 }

 int MyVector::size()const{
	 return m_data.size();
 }
 void MyVector::randWithinRadi(double radius, ProgramMode mode){
	 if(Global::msp_global.get()!=nullptr){
		double r=Global::msp_global->getRandFloat(0,radius,mode);
		randomize(-1,1,mode);
		normalize();
		for(auto&i:m_data) i*=r;
		m_write=true;
	 }else{
		 throw myException("Global::msp_global not initialized@MyVector::randWithinRadi()");
	 }
 }
 void MyVector::randOnRadi(double radius, ProgramMode mode){
	 if(Global::msp_global.get()!=nullptr){
		randomize(-1,1,mode);
		normalize();
		for(auto&i:m_data) i*=radius;
		m_length=radius;
		m_write=false;
	 }else{
		 throw myException("Global::msp_global not initialized@MyVector::randOnRadi()");
	 }
 }

 void MyVector::randOnRadi(double radius, uniform_real_distribution<double> &unif, default_random_engine &gen){

	 randomize(unif, gen, -1, 1);
	 normalize();
	 for (auto&i : m_data) i *= radius;
	 m_length = radius;
	 m_write = false;

 }
 void MyVector::randomize(uniform_real_distribution<double> &unif, default_random_engine &gen, double min, double max){
	 for (auto&i : m_data) i = min+(unif(gen)-unif.min())/(unif.max()-unif.min())*(max - min);
	 m_write = true;
 }
 void MyVector::randomize(double min,double max,ProgramMode mode){	
	 if(Global::msp_global.get()!=nullptr){
		for(auto&i:m_data) i=Global::msp_global->getRandFloat(min,max,mode);	
		m_write=true;
	 }else{
		 throw myException("Global::msp_global not initialized@MyVector::randomize()");
	 }
 }
 void MyVector::norRandomize(ProgramMode mode){
	if(Global::msp_global.get()!=nullptr){
		if(mode==Program_Algorithm)		for(auto&i:m_data) i=Global::msp_global->mp_normalAlg->Next();
		else for(auto&i:m_data) i=Global::msp_global->mp_normalPro->Next();
		m_write=true;
	 }else{
		 throw myException("Global::msp_global not initialized@MyVector::norRandomize()");
	 }
 }
 void MyVector::norRandWithinRadi(double radius, ProgramMode mode){
	if(Global::msp_global.get()!=nullptr){
		double r=0;
		if(mode==Program_Algorithm)		r=fabs(Global::msp_global->mp_normalAlg->NextNonStand(0,radius));
		else r=fabs(Global::msp_global->mp_normalPro->NextNonStand(0,radius));
		randomize(-1,1,mode);
		normalize();
		for(auto&i:m_data) i*=r;
		m_write=true;
	 }else{
		 throw myException("Global::msp_global not initialized@MyVector::randWithinRadi()");
	 }
 }
 MyVector &MyVector::operator =(const vector<double> & v){
	 if(m_data.size()==v.size()){
		 m_data=v;
		 m_write=true;
	 } else{
		 throw myException("size not the same @ MyVector::operator =(vector<double>&)");
	 }
	 return *this;
 }

 MyVector & MyVector::operator *=(double val){
	 for(auto &i:m_data) i*=val;
	 m_write=true;
	 return *this;
 }
 MyVector & MyVector::operator /=(double val){
	for(auto &i:m_data) i/=val;
	  m_write=true;
	 return *this;
 }
 MyVector & MyVector::operator -=(double val){
	for(auto &i:m_data) i-=val;
	m_write=true;
	 return *this;
 }
 MyVector & MyVector::operator +=(double val){
	for(auto &i:m_data) i+=val;
	m_write=true;
	 return *this;
 }
 vector<double>::iterator  MyVector::begin(){
	 return m_data.begin();
 }
 vector<double>::iterator  MyVector::end(){
	 return m_data.end();
 }
 MyVector  MyVector::operator *(double val)const{
	vector<double> v(m_data);
	for(auto&i:v) i*=val;
	return MyVector(v);
 }
 MyVector  MyVector::operator *(double val){
	 vector<double> v(m_data);
	 for (auto&i : v) i *= val;
	 return MyVector(v);
 }
 MyVector  MyVector::operator /(double val){
	vector<double> v(m_data);
	for(auto&i:v) i/=val;
	return MyVector(v);
 }
 MyVector  MyVector::operator -(double val){
 	vector<double> v(m_data);
	for(auto&i:v) i-=val;
	return MyVector(v);
 }
 MyVector  MyVector::operator +(double val){
 	vector<double> v(m_data);
	for(auto&i:v) i+=val;
	return MyVector(v);
 }
 double MyVector::getAngle(MyVector & v){
	 double sum=this->operator*(v);
	 return acos(sum/(length()*v.length()));
 }

void MyVector::push_back(double n){
	m_data.push_back(n);
	m_write=true;
}

double MyVector::getDis(const MyVector&v){
	if(v.size()!=size()) throw myException("size not the same@  MyVector::getDis(const MyVector&v)");
	MyVector vec(m_data);
	vec-=v;
	return vec.length();
}
double MyVector::getDis(const vector<double>&v){
	if(v.size()!=size()) throw myException("size not the same@  MyVector::getDis(const MyVector&v)");
	MyVector vec(m_data);
	vec-=v;
	return vec.length();
}

ifstream &operator>>(ifstream & in,  MyVector&v){
	for(auto &i:v.m_data){
		in>>i;
	}
	v.m_write=true;
	return in;
}
double MyVector::length(){ 
	if(m_write){
		length_();
		m_write=false;
	}
	return m_length;
}
MyVector::MyVector(MyVector&& rhs):m_data(move(rhs.m_data)),m_length(rhs.m_length),m_write(rhs.m_write){
	rhs.m_length = 0;
	rhs.m_write = false;
}

MyVector& MyVector::operator=(MyVector&& rhs){
	m_data = move(rhs.m_data);
	m_length=rhs.m_length;
	m_write=rhs.m_write;
	rhs.m_length = 0;
	rhs.m_write = false;
	return *this;
}

void MyVector::zero(){
	for (auto &i : m_data){
		i=0;
	}
	m_length = 0;
	m_write = false;
}

MyVector MyVector::getPointBetween(const MyVector & v1, double r){
	vector<double> v(m_data);
	for (auto i = 0; i < m_data.size(); ++i){
		v[i] = m_data[i] * r + v1[i] * (1 - r);
	}
	return MyVector(v);
}