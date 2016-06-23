/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 21 September 2011
// Last modified: 12 Dec. 2014
//MyVector.h
#ifndef MYVECTOR_H
#define MYVECTOR_H

#include "definition.h"
#include "TypeVar/typeVar.h"
#include <random>
class Matrix;
class MyVector{								
	vector<double> m_data;
	double m_length;
	bool m_write;
public:
	class MyVectorProxy{
	private: 
		MyVector & m_vec;
		int m_idx;
	public:
		MyVectorProxy(MyVector& v,int idx);
		MyVectorProxy& operator=(const MyVectorProxy& rhs);
		operator double () const;
		MyVectorProxy&  operator=(double val);
		MyVectorProxy& operator+=(double);
	};
	friend class MyVectorProxy;
	MyVector(int s=0);
	MyVector(int s,double val);
	MyVector(int s,double *val);
	MyVector(const vector<double> & v);
	MyVector(const vector<TypeVar> &v);
	MyVector(const MyVector&) = default;
	~MyVector();
	MyVector &operator =(const MyVector & v);
	MyVector &operator =(const vector<double> & v);
	MyVector & operator +=(const MyVector & v);
	MyVector & operator -=(const MyVector & v);
	MyVector & operator *=(double val);
	MyVector & operator /=(double val);
	MyVector & operator -=(double val);
	MyVector & operator +=(double val);
	MyVector  operator *(double val) const;
	MyVector  operator *(double val);
	MyVector  operator /(double val);
	MyVector  operator -(double val);
	MyVector  operator +(double val);
	MyVector getPointBetween(const MyVector & v, double ratio);
	MyVectorProxy operator [](int);
	const double &operator [](int)const;
	double operator *(const MyVector & v);
	void projectionToV(const MyVector &v);
	void normalize();
	int size()const;
	friend class Matrix;
	void randWithinRadi(double radius, ProgramMode mode=Program_Algorithm);
	void randOnRadi(double radius, ProgramMode mode=Program_Algorithm);
	void randOnRadi(double radius, uniform_real_distribution<double> &unif, default_random_engine &gen);
	void randomize(double min=0,double max=1,ProgramMode mode=Program_Algorithm);
	void randomize(uniform_real_distribution<double> &unif, default_random_engine &gen, double min, double max);
	void norRandomize(ProgramMode mode=Program_Algorithm);
	void norRandWithinRadi(double radius, ProgramMode mode=Program_Algorithm);
	vector<double>::iterator  begin();
	vector<double>::iterator  end();
	double getAngle(MyVector & v);
	double length();
	double getDis(const MyVector&v);
	double getDis(const vector<double>&v);
	void push_back(double);
	MyVector(MyVector&& rhs);
	MyVector& operator=(MyVector&& rhs);
	vector<double> & data();
	const vector<double> & data()const;
	void zero();
protected:
	void length_();
	friend ifstream &operator>>(ifstream &,  MyVector&);
};

 inline void MyVector::length_(){
	m_length=0;
	for(auto &i:m_data) m_length+=i*i;
	m_length=sqrt(m_length);
 }
 inline vector<double> & MyVector::data(){
	 return m_data;
 }
 inline const vector<double> & MyVector::data()const{
	 return m_data;
 }
#endif
