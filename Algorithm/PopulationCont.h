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

#ifndef POPULATION_CONT_H
#define POPULATION_CONT_H
#include "Population.h"
#include "../Problem/ContOptima.h"

template<typename>  class MultiPopulationCont;

template<typename ED, typename TypeIndi,typename TypePop=Population<ED,TypeIndi> >
class PopulationCont: public TypePop{
	template<typename> friend class MultiPopulationCont;
protected:
	Solution<ED>  m_center;                        // the center of the population
	double m_initialRadius;						// assumption that a population has a search area, which is defined the average distance of all inidividuals to the pop's center
	double m_curRadius;
	double m_radiusQaulity; /*to evaluate peaks tracked qaulity in terms of population radius*/
	double m_convergThreshold;			// absolute value of radius
	double m_convFactor;				//relative factor 
protected:
	using  TypePop::initialize;
	void initialize(const Solution<ED> & center, double radius, bool mode, bool clearOldBest=true);
	void reInitializeCenter(Optima<ContOptimum> &opt, Solution<ED> &chr);
	double evaluateRadius(double radius);
	ReturnFlag updateBest(const int p);
	ReturnFlag updateBest(const int p,double ratio);
public:
	virtual ~PopulationCont(){ }
	PopulationCont():TypePop(),m_center(),m_initialRadius(0),m_curRadius(0),m_radiusQaulity(0),m_convergThreshold(0.0001),m_convFactor(0.005){
	
	}
	PopulationCont(const int rPopsize,bool mode);
	PopulationCont(const Solution<ED> & center, double radius,const int rPopsize,bool mode);
	PopulationCont( const PopulationCont<ED,TypeIndi,TypePop> &s);
	PopulationCont( Group<ED,TypeIndi> &g);

	virtual void computeCenter();
	virtual void updateCurRadius(bool mode=false);
	virtual void computeInitialRadius();

	void updateRadiusQaulity( int idxPeak, int numPeaksIn,double radius);
	ReturnFlag add(int num,  bool mode,  bool insize);
	void add( TypeIndi *p,bool tranship=true);
	void add(TypePop &s);
	void add(Group<ED,TypeIndi> &g);

	void add( vector<TypeIndi*> &indis);
	void add( vector<unique_ptr<TypeIndi>> &indis);
	void remove(int num,const int *id=0);
	void remove(const vector<int>&);
	int findNearest(int idx,double *rdis,int *sortedIdx,int mode);

	void findNearestPair(int & idx1, int &idx2, int mode);
	void setRadius(double rRadius);
	virtual void increaseDimension();
	virtual void decreaseDimension();
	void updateMemory();
	double GetAvgDistance(bool mode);
	void setConvergThreshold(double value);
	void setCongergFactor(double val);
	bool isConverged();
	bool isConverged(double val);
	Solution<ED> * getNearestBest2Peak(const ED& peak);
	void printToScreen();
	PopulationCont<ED,TypeIndi,TypePop>& operator=(const PopulationCont<ED,TypeIndi,TypePop>& rhs);
	double getInitialRadius();
	Solution<ED>& getCenter();
};

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::initialize(const Solution<ED> & center, double radius, bool mode, bool clearOldBest){
	for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->initialize(center,radius,i,i+1,mode); // index, ID
	this->findWorst();
	this->findBest();
	if(clearOldBest){
		this->m_best.resize(this->m_bestIdx.size());
		for(unsigned i=0;i<this->m_bestIdx.size();i++) this->m_best[i].reset(new Solution<ED>(this->m_pop[this->m_bestIdx[i]]->representative()));
	}else{
		for(unsigned i=0;i<this->m_bestIdx.size();i++){
			this->updateBestArchive(this->m_pop[this->m_bestIdx[i]]->representative());
		}
	}
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::reInitializeCenter(Optima<ContOptimum> &opt, Solution<ED> &chr){
	// get one random point that is not in any areas of opt. NOTICE: no overlap between opts
	Solution<ED> rand;
	rand.initialize();
	for(unsigned int i=0;i<opt.getNumOpt();i++){
		if(rand.getDistance(opt[i].self())<opt[i].m_radius){
			//TODU for non-continious space
			double normal=0;
			double x1=0.,x2=0.,x3=0.;
			for(int d=0;d<GET_NUM_DIM;d++){
				x1=rand[d];x2=opt[i].data()[d];
				normal+=(x1-x2)*(x1-x2);
			}
			normal=sqrt(normal);
				
			for(int d=0;d<GET_NUM_DIM;d++){
				x1=(rand[d]); x2=(opt[i].data()[d]);
				rand[d]=x1+(opt[i].m_radius)*(x1-x2)/normal;
				x3=(rand[d]);
			}
			break;
		}

	}
	chr=rand;
}

template<typename ED, typename TypeIndi,typename TypePop>
double PopulationCont<ED,TypeIndi,TypePop>::evaluateRadius(double radius){
if(Global::msp_global->mp_problem->m_name.find("DYN_CONT_MovingPeak")==string::npos ) return 0;
double minDis;int nearest=-1,numPeaksIn=0;
// get the neaest peak
for(int k=0;k<CAST_PROBLEM_DYN->getNumberofPeak();k++ ){
	if(!CAST_PROBLEM_DYN_CONT->isVisable(k)) continue;
	const double *peak=CAST_PROBLEM_DYN_CONT->getPeak(k);
	vector<double> vdis;
	CodeVReal v(peak,peak+this->m_numDim);
	for(auto &i:this->m_best) vdis.push_back(i->getDistance(v));
	double dis=*min_element(vdis.begin(),vdis.end());

	if(nearest==-1){
		minDis=dis;
		nearest=k;
	}else{
		if(minDis>dis){
			minDis=dis;
			nearest=k;
		}
	}
}
double *p_nearest=const_cast <double *> (CAST_PROBLEM_DYN_CONT->getPeak(nearest));

for(int k=0;k<CAST_PROBLEM_DYN->getNumberofPeak();k++ ){
	if(!CAST_PROBLEM_DYN_CONT->isVisable(k)) continue;
	double *peak=const_cast <double *> (CAST_PROBLEM_DYN_CONT->getPeak(k));
	double dis=0;
	for(int d=0;d<GET_NUM_DIM;d++) dis+=(p_nearest[d]-peak[d])*(p_nearest[d]-peak[d]);

	dis=sqrt(dis);

	if(dis<radius)	numPeaksIn++;

}

updateRadiusQaulity(nearest,numPeaksIn,radius);
return m_radiusQaulity;
}


template<typename ED, typename TypeIndi,typename TypePop>
PopulationCont<ED,TypeIndi,TypePop>::PopulationCont(const int rPopsize,bool mode):TypePop(rPopsize, mode),m_center(),m_radiusQaulity(0),m_convergThreshold(0.0001),m_convFactor(0.005){
	computeCenter();
	computeInitialRadius();	
}

template<typename ED, typename TypeIndi,typename TypePop>
PopulationCont<ED,TypeIndi,TypePop>::PopulationCont(const Solution<ED> & center, double radius,const int rPopsize,bool mode):TypePop(rPopsize),\
	m_center(),m_radiusQaulity(0),m_convergThreshold(0.0001),m_convFactor(0.005){
	initialize(center,radius,mode);
	computeCenter();
	computeInitialRadius();		
}

template<typename ED, typename TypeIndi,typename TypePop>
PopulationCont<ED,TypeIndi,TypePop>::PopulationCont( const PopulationCont<ED,TypeIndi,TypePop> &s):TypePop(s),m_center(s.m_center){	
	m_initialRadius=s.m_initialRadius;
	m_curRadius=s.m_curRadius;
	m_radiusQaulity=s.m_radiusQaulity;
	m_convergThreshold=s.m_convergThreshold;
	m_convFactor=s.m_convFactor;
}

template<typename ED, typename TypeIndi,typename TypePop>
PopulationCont<ED,TypeIndi,TypePop>::PopulationCont( Group<ED,TypeIndi> &g):TypePop(g),m_center(g.getCenter()),m_convergThreshold(0.0001),m_convFactor(0.005){	
	//computeCenter();
	computeInitialRadius();
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::computeCenter(){
	if(this->m_popsize<1) return;
		
	if(Global::msp_global->mp_problem->isProTag(CONT)){
		for( int i=0;i<GET_NUM_DIM;i++){
			double x=0.;
			for(int j=0;j<this->m_popsize;j++){
				Solution<CodeVReal> & chr=this->m_pop[j]->representative();
				x=chr.data()[i]+x;
			}
			m_center.data()[i]=x/this->m_popsize;

		}
			
		m_center.evaluate(false);
		this->updateBestArchive(m_center);
	}

}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::updateCurRadius(bool mode){
		//mode=true for distance between representatice to the center; else between pself to the center
		m_curRadius=0;
		if(this->m_popsize<2) return;

		if(mode){
			for(int j=0;j<this->m_popsize;j++) m_curRadius+=this->m_pop[j]->representative().getDistance(m_center);
		}else{
			for(int j=0;j<this->m_popsize;j++) m_curRadius+=this->m_pop[j]->self().getDistance(m_center);
		}
		m_curRadius=m_curRadius/this->m_popsize;

}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::computeInitialRadius(){
			m_initialRadius=0;
			m_curRadius=0;
			m_radiusQaulity=0;
			if(this->m_popsize<2) return;

			updateCurRadius();
			m_initialRadius=m_curRadius;
			evaluateRadius(m_initialRadius);
	}
		
template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::updateRadiusQaulity( int idxPeak, int numPeaksIn,double radius){
		double aradius= CAST_PROBLEM_DYN_CONT->getAssociateRadius(idxPeak);
		if(radius>aradius){
			m_radiusQaulity=pow(aradius/radius,1+(numPeaksIn-1.)/CAST_PROBLEM_DYN->getNumberofPeak());
		}else{
			m_radiusQaulity=pow(radius/aradius,1+(numPeaksIn-1.)/CAST_PROBLEM_DYN->getNumberofPeak());
		}
	}

template<typename ED, typename TypeIndi,typename TypePop>
ReturnFlag PopulationCont<ED,TypeIndi,TypePop>::add(int num,  bool mode,  bool insize){
	if(num<=0) return Return_Normal;
	//this->m_popsize+=num;	
	ReturnFlag rf=Return_Normal;
	int count=0;
	for(int i=0, k=0;i<num;i++,k++){
		this->m_pop.push_back(move(unique_ptr<TypeIndi>(new TypeIndi())));
		if(insize) rf=this->m_pop[this->m_popsize]->initialize(this->m_best[k%this->m_best.size()]->self(),m_initialRadius,this->m_popsize,this->m_popsize+1,mode);
		rf=this->m_pop[this->m_popsize]->initialize(this->m_popsize,this->m_popsize+1,this->m_popsize+1,mode);
		this->m_popsize++;
		count++;
		if(rf!=Return_Normal) break;
	}
	this->update(count);
	computeCenter();
	updateCurRadius();
	if(m_initialRadius==0) m_initialRadius=m_curRadius;
	return rf;
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::add( TypeIndi *p,bool tranship){
	if(p==nullptr) return;
	TypePop::add(p,tranship);
	computeCenter();
	updateCurRadius();	
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::add(TypePop &s){
	if(s.getPopSize()==0) return;
	TypePop::add(s);
	computeCenter();
	updateCurRadius();	
	//printToScreen();
}

template<typename ED, typename TypeIndi,typename TypePop> 
void PopulationCont<ED,TypeIndi,TypePop>::add(Group<ED,TypeIndi> &g){
	if(g.getSize()==0) return;
	TypePop::add(g);
	computeCenter();
	updateCurRadius();	
}     

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::add( vector<TypeIndi*> &indis){
	if(indis.size()==0) return;
	TypePop::add(indis);
	computeCenter();
	updateCurRadius();	
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::add( vector<unique_ptr<TypeIndi>> &indis){
	if(indis.size()==0) return;
	TypePop::add(indis);
	computeCenter();
	updateCurRadius();	
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::remove(int num,const int *id){
	if(num<=0) return;
	TypePop::remove(num,id);
	if(this->m_popsize==0){
		m_curRadius=0;
		m_initialRadius=0;
		m_radiusQaulity=0;
		return;
	}
	computeCenter();
	updateCurRadius();
	//printToScreen();
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::remove(const vector<int> & id){
	if(id.size()==0) return;
	TypePop::remove(id);
	if(this->m_popsize==0){
		m_curRadius=0;
		m_initialRadius=0;
		m_radiusQaulity=0;
		return;
	}
	computeCenter();
	updateCurRadius();
	//printToScreen();
}

template<typename ED, typename TypeIndi,typename TypePop>
int PopulationCont<ED,TypeIndi,TypePop>::findNearest(int idx,double *rdis,int *sortedIdx,int mode){
	//mode=1: repr to repr; mode=2: self to repre; mode=3 repr to self; mode=4: self to self
    double Min_dis;
	int index;
    int count=0;
	double *dis=new double[this->m_popsize];
	for(int i=0;i<this->m_popsize;i++){
		if(idx==i) {
			dis[i]=0;
			continue;
		}
		switch(mode){
            case 1:
                dis[i]=this->m_pop[idx]->representative().getDistance(this->m_pop[i]->representative());
                break;
            case 2:
                dis[i]=this->m_pop[idx]->self().getDistance(this->m_pop[i]->representative());
                break;
            case 3:
                 dis[i]=this->m_pop[idx]->representative().getDistance(this->m_pop[i]->self());
                 break;
            case 4:
                dis[i]=this->m_pop[idx]->self().getDistance(this->m_pop[i]->self());
                break;
            default:
				throw myException("ensure mode is one of the value in [1,2,3,4]:PopulationCont::findNearest()");
                
		}
	}

	for(int i=0;i<this->m_popsize;i++){
		if(i==idx) continue;
		count++;
        if(count==1) {
            Min_dis=dis[i];
            index=i;
        }else if(dis[i]<Min_dis){
			Min_dis=dis[i];
			index=i;
		}
	}
	if(!sortedIdx){
		if(rdis) *rdis=Min_dis;
		delete [] dis;
		dis=0;
		return index;
	}
	vector<int> idex;
	gQuickSort(dis,this->m_popsize,idex);
	copy(idex.begin(),idex.end(),sortedIdx);
	
	if(rdis)  *rdis=Min_dis;
	delete [] dis;
	dis=0;

	return index;
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::findNearestPair(int & idx1, int &idx2, int mode){

    double mindis,dis;
    for(int i=0;i<this->m_popsize;i++){
       int idx=findNearest(i,&dis,0,mode);
       if(i==0|| mindis>dis){
            idx1=i;
            idx2=idx;
            mindis=dis;
       }
	}
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::setRadius(double rRadius){
	// for algorithms needed to set radius by users
	m_initialRadius=rRadius;
	size_t start, end;
    start=this->m_algPar.str().find("Radius: ");
    for(size_t i=start;i<this->m_algPar.str().size();i++){
        if(this->m_algPar.str()[i]==';') {
            end=i;
            break;
        }
    }
    stringstream ss;
    ss<<"Radius:  "<<m_initialRadius<<";";
	if(start!=string::npos){
		string result=this->m_algPar.str();
		result.replace(start,end-start+1, ss.str());
		 this->m_algPar.str(result);
	}else this->m_algPar<<ss.str();

}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::increaseDimension(){
	TypePop::increaseDimension();
	m_center.increaseDimension();
	computeCenter();
	updateCurRadius();
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::decreaseDimension(){
	TypePop::decreaseDimension();
	m_center.decreaseDimension();
	computeCenter();
	updateCurRadius();
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::updateMemory(){
	if(this->m_popsize==0) return;
	TypePop::updateMemory();
	m_center.evaluate(false);
}

template<typename ED, typename TypeIndi,typename TypePop>
double PopulationCont<ED,TypeIndi,TypePop>::GetAvgDistance(bool mode){
    // mode=true for distance between pself else between representative
    if(this->m_popsize<2) return 0;
    double dis=0;
    for(int i=0;i<this->m_popsize;i++){
        for(int j=i+1;j<this->m_popsize;j++){
            if(mode) dis+=this->m_pop[i]->getDistance(*this->m_pop[j]);
            else dis+=this->m_pop[i]->representative().getDistance(this->m_pop[j]->representative());
        }
    }
    return 2*dis/(this->m_popsize*(this->m_popsize-1));

}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::setConvergThreshold(double value){
	m_convergThreshold=value;
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::setCongergFactor(double val){
	m_convFactor=val;
}

template<typename ED, typename TypeIndi,typename TypePop>
bool PopulationCont<ED,TypeIndi,TypePop>::isConverged(){
	return m_curRadius<=m_convergThreshold;
}

template<typename ED, typename TypeIndi,typename TypePop>
bool PopulationCont<ED,TypeIndi,TypePop>::isConverged(double val){
	return m_curRadius<=val;
}

template<typename ED, typename TypeIndi,typename TypePop>
Solution<ED> * PopulationCont<ED,TypeIndi,TypePop>::getNearestBest2Peak(const ED& peak){
	if(this->m_best.size()==1) return this->m_best[0].get();
	double dis=this->m_best[0]->getDistance(peak);
	unsigned idx=0;
	for(unsigned i=1;i<this->m_best.size();i++){
		double d=this->m_best[i]->getDistance(peak);
		if(d<dis){
			dis=d; idx=i;
		}
	}
	return this->m_best[idx].get();
}

template<typename ED, typename TypeIndi,typename TypePop>
void PopulationCont<ED,TypeIndi,TypePop>::printToScreen(){
	//cout<<this->m_popID<<endl;
	for(int i=0;i<this->m_popsize;i++) this->m_pop[i]->printToScreen();
	cout<<"curdius: "<<m_curRadius<<" Initial radius: "<<m_initialRadius<<" Radius Quality: "<<m_radiusQaulity<<endl;
	cout<<"best: ";
	for(auto&i:this->m_best) i->printToScreen();
	cout<<"center: ";
	m_center.printToScreen();
	cout<<endl;
}

template<typename ED, typename TypeIndi,typename TypePop>
PopulationCont<ED,TypeIndi,TypePop>& PopulationCont<ED,TypeIndi,TypePop>::operator=(const PopulationCont<ED,TypeIndi,TypePop>& rhs){
	if(this==&rhs) return *this;
	TypePop::operator=(rhs);

	m_center=rhs.m_center;                       
	m_initialRadius=rhs.m_initialRadius;						
	m_curRadius=rhs.m_curRadius;
	m_radiusQaulity=rhs.m_radiusQaulity; 
	m_convergThreshold=rhs.m_convergThreshold;			
	m_convFactor=rhs.m_convFactor;	
	return *this;
}
template<typename ED, typename TypeIndi,typename TypePop>
double PopulationCont<ED,TypeIndi,TypePop>::getInitialRadius(){
	return m_initialRadius;
}
template<typename ED, typename TypeIndi,typename TypePop>
Solution<ED> & PopulationCont<ED,TypeIndi,TypePop>::getCenter(){
	return m_center;
}

template<typename ED, typename TypeIndi,typename TypePop>
ReturnFlag PopulationCont<ED,TypeIndi,TypePop>::updateBest(const int p,double ratio){
    Solution<ED> x;
	ReturnFlag r_flag=Return_Normal;
	for(unsigned k=0;k<this->m_best.size();k++){
		for( int j=0;j<this->m_numDim;j++){
			if(ratio<1&&Global::msp_global->mp_uniformAlg->Next()>ratio) continue;
			x=*this->m_best[k];
			x.data()[j]=this->m_pop[p]->data()[j];
			r_flag=x.evaluate();

			if(x>*this->m_best[k])   *this->m_best[k]=x;
			if(r_flag!=Return_Normal) break;
		}
		if(r_flag!=Return_Normal) break;
	}
	return r_flag;

}

template<typename ED, typename TypeIndi,typename TypePop>
ReturnFlag PopulationCont<ED,TypeIndi,TypePop>::updateBest(const int p){
    
	ReturnFlag r_flag=Return_Normal;
	double dis=0;
	vector<double> ratio(this->m_numDim,0);
	for(unsigned k=0;k<this->m_best.size();k++){
		dis= this->m_best[k]->getDistance(*this->m_pop[p],DIS_MANHATTAN);
		if(dis>0){
			Solution<ED> & self= *this->m_pop[p];
			Solution<ED> & best=*this->m_best[k];
			for( int j=0;j<this->m_numDim;j++){
				ratio[j]=1-fabs((double)(self.data()[j]-best.data()[j]))/dis;
			}
		}else{
			continue;
		}

		Solution<ED> x;
		for( int j=0;j<this->m_numDim;j++){
			double r=Global::msp_global->mp_uniformAlg->Next();
			if(r>ratio[j]) continue;
			x=*this->m_best[k];
			Solution<ED> & self= *this->m_pop[p];
			x.data()[j]=self.data()[j];
			r_flag=x.evaluate();
			if(x>*this->m_best[k])   *this->m_best[k]=x;
			if(r_flag!=Return_Normal) break;
		}
		if(r_flag!=Return_Normal) break;
	}
	return r_flag;
}


#endif