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

#ifndef Matrix_H
#define Matrix_H

#include "myVector.h"

#ifdef use_namespace
namespace Matrix{
#endif

class Matrix{								// *****************Orthogonal rotation matrix***********************
private:
	int m_col,m_row;									// matrix size
	vector<MyVector> m_vec;
public:
	Matrix(int dim=0);		
	Matrix(const int c,const int r);
	~Matrix();

	Matrix & operator*=(const Matrix & m);
	Matrix & operator*=(const double x);
	Matrix & operator=(const Matrix & m);
	bool identity();
	bool isIdentity();
	void setRotationAngle(const int r,const int c,const double angle);
	void generateRotationMatrix(const double CondiNum,ProgramMode rMode);
	void randomize(ProgramMode rMode);
	void orthonormalize();
	void setToZero();
	void setDataRow(const double *d, const int c,const int r=1);
	void setDataCol(const double *d, const int r,const int c=1);
	void setData(const double * const * d);
	void diagonalize(const double CondiNum );
	void transpose();
	void inverse();
	MyVector & operator[](const int idx);
	const MyVector & operator[](const int idx)const;
	void Read_Data(ifstream &in){
		for(int i=0;i<m_row;i++){
			in>>m_vec[i];
		}
	};
	// used for debug
	void Print(ofstream & out){
		for(int i=0;i<m_row;i++){
			for(int j=0;j<m_col;j++)
				out<<m_vec[i][j]<<" ";
			out<<endl;
		}
		//out<<endl;
	};
	void print();
	void resize(const int row, const int col);
	
private:
        void freeMemory();
        void allocateMemory(const int r, const int c);
};
#ifdef use_namespace
}
#endif

#endif
