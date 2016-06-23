#include "FAMFPopDE.h"



ReturnFlag FAMFPopDE::evolve(){
    if(this->m_popsize<5){
		throw myException("the population size cannot be smaller than 5@DEPopulation<TypeDEIndi>::evolve()");
    }
	ReturnFlag r_flag=Return_Normal;
    for(int i=0;i<this->m_popsize;i++){
        mutate(i);
        this->m_pop[i]->recombine(m_CR);
    }

    this->updateIDnIndex();
    for(int i=0;i<this->m_popsize;i++){
        r_flag=this->m_pop[i]->select();
		this->updateBestArchive(this->m_pop[i]->self());
		if(r_flag!=Return_Normal) break;

    }
	if(r_flag==Return_Normal) {
		computeCenter();
		this->updateCurRadius(true);
		this->m_iter++;
	}
    return r_flag;
}