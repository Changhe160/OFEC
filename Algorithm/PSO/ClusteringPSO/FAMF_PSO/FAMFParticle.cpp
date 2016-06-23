
#include "FAMFParticle.h"
#include "../../../../Utility/myVector.h"
#include "../../../Clustering/FAMFDera.h"

FAMFParticle::FAMFParticle():Particle(){
		
}
FAMFParticle::FAMFParticle( const Solution<CodeVReal> &chr):Particle(chr){
		
}

void FAMFParticle::initializeVelocityAftClustering(){
	for(int i=0;i<GET_NUM_DIM;i++){
		m_vel[i]=(Global::msp_global->mp_uniformAlg->Next()-0.5)*(m_vMax[i].m_max-m_vMax[i].m_min);
	}
}
// Brownian movements to dealwith noisy environments for the gbest particle
ReturnFlag FAMFParticle::brwonianMove(double radius){
	//self().printToScreen();
	for(int i=0;i<GET_NUM_DIM;i++){
		data()[i]+=Global::msp_global->mp_normalAlg->NextNonStand(0,radius);
	}
	self().validate();
	ReturnFlag rf=Return_Normal;
	if(FAMFDerating::derateFitness(self())==false)	rf=self().evaluate();
	
	return rf;
}

ReturnFlag FAMFParticle::cauchyMove(double radius){
	//self().printToScreen();
	for(int i=0;i<GET_NUM_DIM;i++){
		if(radius<0){
			double l,u;
			CAST_PROBLEM_CONT->getSearchRange(l,u,i);
			data()[i]+=Global::msp_global->mp_cauchyAlg->NextNonStand(0,(u-l)/2);
			}else{
			data()[i]+=Global::msp_global->mp_cauchyAlg->NextNonStand(0,radius);				
		}
	}
	self().validate();
	ReturnFlag rf=Return_Normal;
	if(FAMFDerating::derateFitness(self())==false)	rf=self().evaluate();
	
	
	return rf;
}
ReturnFlag FAMFParticle::move( const Solution<CodeVReal> &lbest,double w, double c1, double c2){
	double u,l,x;

	for( int j=0;j<GET_NUM_DIM;j++){
		CAST_PROBLEM_CONT->getSearchRange(l,u,j);
		double r1=Global::msp_global->mp_uniformAlg->Next();
		double r2=Global::msp_global->mp_uniformAlg->Next();
		x=data()[j];
		m_vel[j]=w*m_vel[j]+c1*r1*(m_pbest.data()[j]-x)+	c2*r2*((lbest.data()[j])-(x));
		
		if(m_vel[j]>m_vMax[j].m_max)	m_vel[j]=m_vMax[j].m_max;
		else if(m_vel[j]<m_vMax[j].m_min)		m_vel[j]=m_vMax[j].m_min;

		data()[j]=x+m_vel[j];
	}
	self().validate();
	ReturnFlag rf=Return_Normal;
	if(FAMFDerating::derateFitness(self())==false)	rf=self().evaluate();
	
	return rf;
}
