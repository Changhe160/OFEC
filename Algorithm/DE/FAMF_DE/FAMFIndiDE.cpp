
#include "FAMFIndiDE.h"
#include "../../../Global/global.h"
#include "../../Clustering/FAMFDera.h"
#include "../../../Problem/ContinuousProblem.h"
FAMFIndiDE::FAMFIndiDE():DEIndividual(){
		
}
FAMFIndiDE::FAMFIndiDE( const Solution<CodeVReal> &chr):DEIndividual(){
	self()=chr;
}


// Brownian movements to dealwith noisy environments for the gbest particle
ReturnFlag FAMFIndiDE::brwonianMove(double radius){
	for(int i=0;i<GET_NUM_DIM;i++){
		data()[i]+=Global::msp_global->mp_normalAlg->NextNonStand(0,radius);
	}
	self().validate();
	ReturnFlag rf=Return_Normal;
	if(FAMFDerating::derateFitness(self())==false)	rf=self().evaluate();
	
	return rf;
}
ReturnFlag FAMFIndiDE::cauchyMove(double radius){
	//self().printToScreen();
	for(int i=0;i<GET_NUM_DIM;i++){		
		if(radius<=0){
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
ReturnFlag FAMFIndiDE::select(){
	ReturnFlag rf=Return_Normal;
	if(FAMFDerating::derateFitness(m_pu)==false)	rf=m_pu.evaluate();

    if(m_pu>self()) self()=m_pu;
	return rf;
}
