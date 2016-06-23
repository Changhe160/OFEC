/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 19 Jan 2015
// Last modified:
#ifndef PROGRESS_H
#define PROGRESS_H

#include "../../../Global/global.h"
class Progress
{
    public:
        int m_numSelected;
        int m_numSuccess;
        double m_rewards;
        double m_ratio;
        double m_minRatio;
    public:
        Progress() {}
        virtual ~Progress() {}
	void initialize( const int num_ls_operators ){
			m_ratio=1./(num_ls_operators);
			m_minRatio=0.001;
			m_numSelected=0;
			m_numSuccess=0;
			m_rewards=0;
	}

	static void updateProgress(Global *glob, vector<Progress> &prog,const int num_ls_operators, double p_factor=0.9){
            //update selection probability of all mutation operators
                double sum1=0,sum2=0;
				vector<double> temp_value(num_ls_operators);

                for(int i=0;i<num_ls_operators;i++) sum1+=prog[i].m_rewards;

                double max_ratio,m_minRatio;
                max_ratio=m_minRatio=prog[0].m_ratio;
                int index_max=0,index_min=0;
                for(int i=1;i<num_ls_operators;i++){
                    if(max_ratio<prog[i].m_ratio) {
                        max_ratio=prog[i].m_ratio;
                        index_max=i;
                    }
                    if(m_minRatio>prog[i].m_ratio) {
                        m_minRatio=prog[i].m_ratio;
                        index_min=i;
                    }
                }
                if(max_ratio!=m_minRatio&& prog[index_max].m_numSuccess==0)	prog[index_max].m_ratio*=p_factor;

                for(int i=0;i<num_ls_operators;i++){
                    double	alpha=glob->mp_uniformAlg->Next();
                    double t=prog[i].m_numSelected>0?(double)prog[i].m_numSuccess/prog[i].m_numSelected:0;
                    if(sum1>0)
                        temp_value[i]=alpha*t+(1-alpha)* prog[i].m_rewards/sum1+prog[i].m_ratio;
                    else
                        temp_value[i]=prog[i].m_ratio;
                    sum2+=temp_value[i];
                }


                for(int i=0;i<num_ls_operators;i++)
                    prog[i].m_ratio=temp_value[i]*(1-num_ls_operators*prog[i].m_minRatio)/sum2+prog[i].m_minRatio;

    }

    static unsigned int getAction(Global *glob,const unsigned int num_actions,const vector<Progress> &p){
               double sum, pick;
               unsigned int i;
                pick = glob->mp_uniformAlg->Next();
                sum = 0;
                for(i = 0; (sum <= pick) && (i < num_actions); i++)
                    sum += p[i].m_ratio;

                return(i-1);
    }

};

#endif // PROGRESS_H
