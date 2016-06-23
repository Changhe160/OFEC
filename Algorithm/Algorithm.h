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
#ifndef ALGORITHM_H
#define ALGORITHM_H
#include "../Utility/definition.h"
#include "../Utility/TypeVar/typeVar.h"
#include "Termination.h"
class Algorithm{
    protected:
       const int m_algID;
	   int m_numDim;
	   unique_ptr<Termination> m_term;
    public:
        string m_name;
        stringstream m_algPar;
        Algorithm(const int rID, string rName);
        virtual ~Algorithm()=0;
		virtual ReturnFlag run_();
		ReturnFlag run();
        int getID(){return m_algID;}
		Algorithm & operator=(const Algorithm& rhs);
		virtual bool ifTerminating();
		bool ifTerminated();
protected:
	virtual ReturnFlag evolve(){return Return_Normal;};
};

#endif // ALGORITHM_H
