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
// Last modified: 21 Nov 2014
#include "SPSO11.h"
#include "../../../Global/global.h"

#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif
extern mutex g_mutexStream;

void ParticleSPSO11::initializeVelocity() {
	double u, l;
	for (int i = 0; i<GET_NUM_DIM; i++) {
		CAST_PROBLEM_CONT->getSearchRange(l, u, i);
		m_vel[i] = Global::msp_global->getRandFloat(l - m_x[i], u - m_x[i]);
	}
}

ReturnFlag ParticleSPSO11::move(double w, double c1, double c2, const Solution<CodeVReal> *lbest,bool clamping){
    double u,l;
	vector<double> P(GET_NUM_DIM),L(GET_NUM_DIM),C(GET_NUM_DIM);

	for (int j = 0; j < GET_NUM_DIM; j++) {
		P[j] = m_x[j] + c1*Global::msp_global->mp_uniformAlg->Next()*((m_pbest.data()[j]) - m_x[j]);
		if (lbest) {
			L[j] = m_x[j] + c2*Global::msp_global->mp_uniformAlg->Next()*(lbest->data()[j] - m_x[j]);
			C[j] = (P[j] + L[j] + m_x[j]) / 3;
		}
		else {
			C[j] = (P[j] + m_x[j]) / 2;
		}
		
	}


	double dis = 0;
	for (int j = 0; j < GET_NUM_DIM; j++) dis += (C[j] - m_x[j])*(C[j] - m_x[j]);
	dis = sqrt(dis);

	MyVector v(GET_NUM_DIM);
	v.randWithinRadi(dis);
	v += C;

	for( int j=0;j<GET_NUM_DIM;j++){
		CAST_PROBLEM_CONT->getSearchRange(l,u,j);
		double x= w*m_vel[j] + v[j];
		m_vel[j] = x - m_x[j];
		m_x[j]=x;

		if(clamping){
			if(m_x[j]>u){
				m_x[j]=u;
				m_vel[j]*=-0.5;
			}else if(m_x[j]<l){
				m_x[j]=l;
				m_vel[j]*=-0.5;
			}
		}	
	}
	
	if(clamping) return self().evaluate();
	
	self().validate();
	return self().evaluate();
}


SPSO11::SPSO11():Swarm<CodeVReal,ParticleSPSO11>()
{
    //ctor
}

SPSO11::SPSO11(ParamMap& v):Swarm<CodeVReal,ParticleSPSO11>(v[param_popSize],v[param_evalCountFlag]),m_impr(0){

	m_name="SPSO11";
    m_algPar<<"Population size: "<<m_popsize;
	initializePara(1. / (2 * log(2)), 0.5 + log(2), 0.5 + log(2));
	m_p=  1 - pow(1 - 1. / (m_popsize), 3);
}

ReturnFlag SPSO11::run_(){
	ReturnFlag rf=Return_Normal;

	#ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
	if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	while(!ifTerminating()){
		/*g_mutexStream.lock();
		cout<<Global::msp_global->m_runId<<" "<<Global::msp_global->mp_problem->getEvaluations()<<" "<<m_best[0]->obj(0)<<endl;
		g_mutexStream.unlock();*/
		#ifdef OFEC_DEMON
			vector<Algorithm*> vp;	
			vp.push_back(this);	
			msp_buffer->updateBuffer_(&vp);
		#endif
        rf=evolve();

		if(rf!=Return_Normal) handleReturnFlag(rf);

		#ifdef OFEC_CONSOLE
		if(mMultiModal::getPopInfor()){
			int peaksf;		
			peaksf=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
			mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
				Global::msp_global->m_totalNumIndis,1,peaksf,\
				CAST_PROBLEM_CONT->getGOpt().getNumOpt(),0,0,0,0,\
				0,0,CAST_PROBLEM_CONT->getGOpt().isAllFound());
		}
		#endif

		if(rf==Return_Terminate) break;
    }
	 
     return rf;
}
ReturnFlag SPSO11::evolve(){

    //lbest model with an adaptive random topology
	
	ReturnFlag rf=Return_Normal;
	if(m_popsize<1) return rf;
	
	vector<int> rindex(m_popsize);
	Global::msp_global->initializeRandomArray(rindex,m_popsize);	//generate a permutation of particle index

	setNeibourhood();

	m_impr=0;
	for(int i=0;i<m_popsize;i++){
		Solution<CodeVReal> *l = &neighborBest(rindex[i]);

		if (l != &m_pop[rindex[i]]->m_pbest)
			rf = m_pop[rindex[i]]->move(m_W, m_C1, m_C2, l, true);
		else rf = m_pop[rindex[i]]->move(m_W, m_C1, m_C2, 0, true);

		
		if(m_pop[rindex[i]]->self()>(m_pop[rindex[i]]->m_pbest)){
			m_pop[rindex[i]]->m_pbest=m_pop[rindex[i]]->self();
			if(updateBestArchive(m_pop[rindex[i]]->self())){
				m_impr++;
			}
		}

		handleReturnFlag(rf);
		HANDLE_RETURN_FLAG(rf)
      
	}
	updateCurRadius();
	m_iter++;
	return rf;
}


void SPSO11::setNeibourhood() {
	for (int i = 0; i<m_popsize; i++) {
		if (m_impr == 0) {
			for (int k = 0; k<m_popsize; k++) {
				if (Global::msp_global->mp_uniformAlg->Next()<m_p) m_link[i][k] = true; //probabilistic method
				else
					m_link[i][k] = false;
			}
			m_link[i][i] = true;
		}
	}
}


Solution<CodeVReal>& SPSO11::neighborBest(int idx) {
	int l = idx;
	for (int k = 0; k<m_popsize; k++) {
		if (m_link[idx][k] == true && m_pop[k]->m_pbest>(m_pop[l]->m_pbest)) {
			l = k;
		}
	}
	return m_pop[l]->m_pbest;
}