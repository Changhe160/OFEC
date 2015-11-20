#include "CPSOSubSwarm.h"

CPSOParticle::~CPSOParticle(){

}

CPSOParticle::CPSOParticle():Particle(){

}
CPSOParticle::CPSOParticle( const Solution<CodeVReal> &chr):Particle(chr){
		
}

///**** Begin of definition of member functions of class CPSOSubSwarm ***///
CPSOSubSwarm::CPSOSubSwarm(int size, bool flag):Swarm<CodeVReal,CPSOParticle>(size,flag){

}
CPSOSubSwarm::CPSOSubSwarm(Group<CodeVReal,CPSOParticle> &g):Swarm<CodeVReal,CPSOParticle>(g){

}

ReturnFlag CPSOSubSwarm::evolve(){
	if(this->m_popsize<1) return Return_Normal;
	Solution<CodeVReal> t;
	ReturnFlag r_flag=Return_Normal;

	for(int i=0;i<this->m_popsize;i++){
		t=this->m_pop[i]->self();
		 bool flag=false;  
        r_flag=this->m_pop[i]->moveBound(this->m_pop[i]->m_pbest,this->getNearestBest(this->m_pop[i]->self()),m_W,m_C1,m_C2);//
			
		if(this->m_pop[i]->self()>this->m_pop[i]->m_pbest){
			this->m_pop[i]->m_pbest=this->m_pop[i]->self();
			flag=this->updateBestArchive(this->m_pop[i]->self());
		}

		if(r_flag!=Return_Normal) break;
		if(!flag&&this->m_pop[i]->self()>t){
			r_flag=updateBest(i,1.0);
		}
		if(r_flag!=Return_Normal)  {  break;}
		
	}
	if(r_flag==Return_Normal){
		this->m_evoNum++;
		double ratio=Global::msp_global->mp_problem->getEvaluations()%CAST_PROBLEM_DYN->getChangeFre()/static_cast<double>(CAST_PROBLEM_DYN->getChangeFre());
		m_W=m_maxW-(m_maxW-m_minW)*ratio;	
	}
	computeCenter();
	updateCurRadius(true);

	return r_flag;
}


