#include "FAMFPopPSO.h"


FAMFPopPSO::FAMFPopPSO( Group<CodeVReal,FAMFParticle> &gr):FAMFPop<FAMFParticle,Swarm<CodeVReal,FAMFParticle>>(gr){
	initializePara(0.7298f,1.496f,1.496f,true);
	initilizeVelAftClustering();

}

void FAMFPopPSO::initilizeVelAftClustering(){
	for(int i=0;i<this->m_popsize;i++){
		m_pop[i]->setVmax(-2.5*m_initialRadius,2.5*m_initialRadius);
		m_pop[i]->initializeVelocityAftClustering();
	}
}

ReturnFlag FAMFPopPSO::evolve(){
	if(this->m_popsize<1||m_flag[0]) return Return_Normal;
	ReturnFlag r_flag=Return_Normal;


	for(int i=0;i<this->m_popsize;i++){	
		// using gbest mode
		//if(Global::msp_global->mp_problem->getEvaluations()>=9175) 
		//	getchar();
		//m_pop[i]->printToScreen();
		r_flag=this->m_pop[i]->move(this->m_pop[i]->m_pbest,this->getNearestBest(this->m_pop[i]->self()),m_W,m_C1,m_C2);//
		//m_pop[i]->printToScreen();
		if(this->m_pop[i]->self()>this->m_pop[i]->m_pbest){
			this->m_pop[i]->m_pbest=this->m_pop[i]->self();
			this->updateBestArchive(this->m_pop[i]->self());
		}
		if(r_flag!=Return_Normal) break;
	}
	if(r_flag==Return_Normal){
		this->computeCenter();		
		this->updateCurRadius(true);
		this->m_evoNum++;
	}

	return r_flag;
}

