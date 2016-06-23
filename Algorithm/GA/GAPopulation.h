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
// Created: 07 Nov 2011
// Last modified:

#ifndef GAPOPULATION_H
#define GAPOPULATION_H

#include "../Population.h"

template<class ED, class GAIndi>
class GAPopulation : public Population<ED,GAIndi>
{
    public:
        GAPopulation<ED,GAIndi>& operator=(const GAPopulation<ED,GAIndi>& rhs);
        GAPopulation();
        virtual ~GAPopulation();

        GAPopulation(int rPopsize,bool mode=true);
        GAPopulation(const GAPopulation<ED,GAIndi> &rhs);
        GAPopulation(Group<ED,GAIndi> &g);
        
        void setDefaultPar();
        void add( GAIndi &p);
        void add(Population<ED,GAIndi> &s);
        ReturnFlag add( int num,  bool mode=true);

		void add(Group<ED,GAIndi> &g);
		void add( vector<GAIndi*> &indis);
		void add( vector<unique_ptr<GAIndi>> &indis);
		void add( GAIndi *indis,bool tranship=true);

		void remove( int num,const int *id=0);
		void remove(const vector<int> &id);

        void setTournSize(const int size);
        void setSelStrategy(GASelectionStrategy sel);
        void setMutationProbability(const double rw);
		void setMutationStrategy(GAMutationStrategy str);

		void setXoverStrategy(GAXoverStrategy str);
        void setMutationStep(const double step);
        void setXoverProbability(const double rw);
protected:
		ReturnFlag evolve();
	    void mapObj2Fitness();
        void rouletteWheel( );
        void makeShuffle( int *shufflearray, const int n );
        // set things up for tournament selection without replacement
        void preTselectWoRep( int *shufflearray, int &pick ) ;
        // tournment selection without replacement
        int getTWinner( int *shuffle, int &pick );
        // tournament selection without replacement
        // the individuals that survive reproduction are placed in the mating pool
        void tSelectWoRep(  );
    protected:
        double m_sumFit;
        int m_tournSize;
        GASelectionStrategy m_selStr;    // selection strategy
        vector<int> m_pool;       // mating pool
};

template<class ED, class GAIndi>
GAPopulation<ED,GAIndi>& GAPopulation<ED,GAIndi>::operator=(const GAPopulation<ED,GAIndi>& rhs){
    if(this==&rhs) return *this;
    Population<ED,GAIndi>::operator=(rhs);
    m_sumFit=rhs.m_sumFit;
    m_selStr=rhs.m_selStr;
    m_tournSize=rhs.m_tournSize;
	m_pool=rhs.m_pool;
    return *this;
}
template<class ED, class GAIndi>
GAPopulation<ED,GAIndi>::GAPopulation():Population<ED,GAIndi>(){
    setDefaultPar();
}
template<class ED, class GAIndi>
GAPopulation<ED,GAIndi>::~GAPopulation(){
    m_pool.clear();
}
template<class ED, class GAIndi>
GAPopulation<ED,GAIndi>::GAPopulation(int rPopsize,bool mode):Population<ED,GAIndi>(rPopsize,mode){
	m_pool.resize(this->m_popsize);
    setDefaultPar();
    mapObj2Fitness();

}
template<class ED, class GAIndi>
GAPopulation<ED,GAIndi>::GAPopulation(const GAPopulation<ED,GAIndi> &rhs):Population<ED,GAIndi>(rhs){
    m_sumFit=rhs.m_sumFit;
    m_selStr=rhs.m_selStr;
    m_tournSize=rhs.m_tournSize;
	m_pool=rhs.m_pool;
}
template<class ED, class GAIndi> 
GAPopulation<ED,GAIndi>::GAPopulation(Group<ED,GAIndi> &g):Population<ED,GAIndi>(g){
	m_pool.resize(this->m_popsize);
    setDefaultPar();
    mapObj2Fitness();
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setDefaultPar(){
    setSelStrategy(SEL_TOURNAMENT);
    setTournSize(3);
          
	setMutationProbability(1./GET_NUM_DIM);
    setXoverProbability(0.6);
	setMutationStrategy(MUTAT_POLYNOMIAL);

}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::mapObj2Fitness(){
    // map objective function values to fitness values.
    // see Goldberg's book, page 75-76
    double maxboj=this->getMaxObj(0);
    double minobj=this->getMinObj(0);
    m_sumFit=0;
        for(int i=0;i<this->m_popsize;i++){
			if(Global::msp_global->mp_problem->getOptType()==MIN_OPT){
            this->m_pop[i]->m_fitness=maxboj-this->m_pop[i]->self().obj(0);
        }else{
            this->m_pop[i]->m_fitness=this->m_pop[i]->self().obj(0)+fabs(minobj);
        }
        m_sumFit+=this->m_pop[i]->m_fitness;
        }
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::rouletteWheel( ){
//---------------------------------------------------
// roulette wheel selection
// see Goldberg, "Genetic Algorithms in Search, Optimization,
// and Machine Learning", 1989
// the individuals that survive reproduction are placed in the mating pool

    double rndpoint, partsum;  // random point on wheel, partial sum
    int    j;                  // population index

    for( int i=0; i< this->m_popsize; i++ ){
        partsum = 0; j = -1;   // Zero out counter and accumulator
        // Wheel point calc. uses random number [0,1]
        rndpoint = Global::msp_global->mp_uniformAlg->Next() * m_sumFit;
        // find wheel slot
        do {
        j++;
        partsum += this->m_pop[i]->m_fitness;
        } while( (partsum < rndpoint) && (j< this->m_popsize -1 ));
        m_pool[i] = j;
    }
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::makeShuffle( int *shufflearray, const int n ){
// tournament selection without replacement
// make a random n-permutation of the numbers 0,1,2,...,n-1
    int i;

    // initialize
    for( i=0; i<n; i++ ) shufflearray[i] = i;
    // shuffle
    for( i=0; i<n-1; i++ ){
		int other = Global::msp_global->getRandInt( i, n-1 );
        // swap( shufflearray[other], shufflearray[i] );
        int temp = shufflearray[other];
        shufflearray[other] = shufflearray[i];
        shufflearray[i] = temp;
    }
}

// set things up for tournament selection without replacement
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::preTselectWoRep( int *shufflearray, int &pick ) {
    pick = 0;
    makeShuffle( shufflearray, this->m_popsize );
}

// tournment selection without replacement
template<class ED, class GAIndi>
int GAPopulation<ED,GAIndi>::getTWinner( int *shuffle, int &pick ){

    if (pick+m_tournSize > this->m_popsize) preTselectWoRep( shuffle, pick );
    int winner = shuffle[pick];
    for( int i=pick+1; i< pick+m_tournSize; i++ ){
    if (this->m_pop[ shuffle[i] ]->m_fitness >= this->m_pop[ winner ]->m_fitness)
		winner = shuffle[i];
	}
    pick += m_tournSize;
    return winner;
}

// tournament selection without replacement
// the individuals that survive reproduction are placed in the mating pool
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::tSelectWoRep(  ){
	if(this->m_popsize<m_tournSize){
		for( int i=0; i< this->m_popsize; i++ ) 		m_pool[i] =i;
		return;
	}
    int *shufflearray;
    shufflearray = new int[this->m_popsize];
    int pick;
    preTselectWoRep( shufflearray, pick );
    for( int i=0; i< this->m_popsize; i++ )           {
        m_pool[i] =getTWinner( shufflearray, pick );
    }
    delete [] shufflearray;
	shufflearray=0;
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::add( GAIndi &p){
    Population<ED,GAIndi>::add(p);
    m_pool.resize(this->m_popsize);
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::add(Population<ED,GAIndi> &s){
    Population<ED,GAIndi>::add(s);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi> 
void GAPopulation<ED,GAIndi>::add(Group<ED,GAIndi> &g){
    Population<ED,GAIndi>::add(g);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
ReturnFlag GAPopulation<ED,GAIndi>::add( int num,  bool mode){
	ReturnFlag rf=  Population<ED,GAIndi>::add(num,mode);
    m_pool.resize(this->m_popsize);
	return rf;
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::add( vector<GAIndi*> &indis){
	Population<ED,GAIndi>::add(indis);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::add( GAIndi *indis,bool tranship){
	Population<ED,GAIndi>::add(indis,tranship);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::add( vector<unique_ptr<GAIndi>> &indis){
	Population<ED,GAIndi>::add(indis);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::remove( int num,const int *id){
    Population<ED,GAIndi>::remove(num,id);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::remove( const vector<int> & id){
    Population<ED,GAIndi>::remove(id);
    m_pool.resize(this->m_popsize);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setTournSize(const int size){
    m_tournSize=size;
    if(m_selStr==SEL_TOURNAMENT){
        size_t start, end;
        start=this->m_algPar.str().find("tournament size:");
        for(size_t i=start;i<this->m_algPar.str().size();i++){
            if(this->m_algPar.str()[i]==';') {
                end=i;
                break;
            }
        }
        stringstream ss;
        ss<<"tournament size:"<<m_tournSize<<"; ";
        if(start!=string::npos){
            string result=this->m_algPar.str();
            result.replace(start,end-start+1, ss.str());
                this->m_algPar.str(result);
        }else this->m_algPar<<ss.str();
    }
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setSelStrategy(GASelectionStrategy sel){
    m_selStr=sel;
    size_t start, end;
    start=this->m_algPar.str().find("Selection strategy:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    switch(m_selStr){
        case SEL_TOURNAMENT:
        ss<<"Selection strategy:tournament selection; ";
        break;
        case SEL_ROULETTE_WHEEL:
        ss<<"Selection strategy:roulette wheel selection; ";
        break;

    }
    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setMutationProbability(const double rw){

    size_t start, end;
    start=this->m_algPar.str().find("Mutation probability:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Mutation probability:"<<rw<<"; ";
    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();

    for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->setMutationProbability(rw);
}
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setMutationStrategy(GAMutationStrategy str){

    size_t start, end;
    start=this->m_algPar.str().find("Mutation strategy:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
	switch(str){
        case MUTAT_POLYNOMIAL:
        ss<<"Mutation strategy:polynomial distribution; ";
        break;
        case MUTAT_NORMAL:
        ss<<"Mutation strategy:normal distribution; ";
        break;
		case MUTAT_COMB:
		ss<<"Mutation strategy:random combinatorial; ";
        break;
    }

    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();

    for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->setMutationStrategy(str);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setXoverStrategy(GAXoverStrategy str){

    size_t start, end;
    start=this->m_algPar.str().find("Xover strategy:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
	switch(str){
        case XOVER_ARITHMETIC:
        ss<<"Xover strategy:arithmetic; ";
        break;
		case XOVER_SINGLEPOINT:
        ss<<"Xover strategy:single point; ";
        break;

    }

    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();

    for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->setXoverStrategy(str);
}

template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setMutationStep(const double step){
        size_t start, end;
    start=this->m_algPar.str().find("Mutation step:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Mutation step:"<<step<<"; ";
    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();

    for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->setMutationStep(step);
}
        
template<class ED, class GAIndi>
void GAPopulation<ED,GAIndi>::setXoverProbability(const double rw){

    size_t start, end;
    start=this->m_algPar.str().find("Crossover probability:");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Crossover probability:"<<rw<<"; ";
    if(start!=string::npos){
        string result=this->m_algPar.str();
        result.replace(start,end-start+1, ss.str());
            this->m_algPar.str(result);
    }else this->m_algPar<<ss.str();

    for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->setXoverProbability(rw);
}

template<class ED, class GAIndi>
ReturnFlag GAPopulation<ED,GAIndi>::evolve(){
    int i;
    //  selection
    switch(m_selStr){
        case SEL_ROULETTE_WHEEL:
        rouletteWheel();
        break;
        case SEL_TOURNAMENT:
        tSelectWoRep();
        break;
    }

	for( i=0; i< this->m_popsize; i++ )            this->m_pop[i]->m_backup=this->m_pop[i]->self();

    for( i=0; i< this->m_popsize; i++ ){
		this->m_pop[i]->self()=this->m_pop[m_pool[i]]->m_backup;
    }

    //
    //  crossover
    //
    i=0;
    while (i<this->m_popsize){
        this->m_pop[i]->xover(*this->m_pop[i+1]);
        i+=2; // next pair of mates
		}
    //
    //  mutation
    //
    for( i=0; i< this->m_popsize; i++ ) this->m_pop[i]->mutate();
    //
    //  compute objective function value of offspring
    //
	ReturnFlag rf=Return_Normal;
    for( i=0; i< this->m_popsize; i++ ){
		rf=this->m_pop[i]->self().evaluate();
		if(rf!=Return_Normal) return rf;
    }
	this->findBest();
	for(unsigned i=0;i<this->m_bestIdx.size();i++){
		this->updateBestArchive(this->m_pop[this->m_bestIdx[i]]->representative());
	} 
    mapObj2Fitness();
	return rf;
}
#endif // GAPOPULATION_H
