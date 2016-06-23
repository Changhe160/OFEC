#ifndef FAMFRAMEWORK_H
#define FAMFRAMEWORK_H

#include "FAMFPop.h"
#include "../MultiPopulationCont.h"

#ifdef OFEC_DEMON
#include "../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

#ifdef OFEC_DEBUG_
static  vector<vector<int>> avgIndi(30);
#endif

template<typename TypeIndi,typename TypeMainPop,typename TypeSubPop>
class FAMF: public TypeMainPop, public MultiPopulationCont<TypeSubPop> {
private:		
	int mc_maxIndis,mc_minIndis,mc_stepIndis;
	//maximum, minimum number of indis allowed in the search space;mc_stepIndis: number of indis to be added or removed when necessary;
	//between two successive diversity increase operation
	const int mc_offPeak; //threshold of the difference of pop numbers between two successive diversity increase
	int m_preIndis,m_nextIndis,m_prePops;
	// total individual at prevous and next diversity increase
	double ms_maxPops;    //largest number of pops since last diversity increase
	vector<unique_ptr<TypeIndi>> mv_convered;
	unique_ptr<TypeMainPop>  mp_single;
	vector< vector<double> > mv_indis;
	Cluster<CodeVReal,TypeIndi> m_clst;
	double m_hibernatingRadius;
	int m_minNumIndis;
	Optima<CodeVReal,ContOptimum> m_optFound;
public:
	FAMF(ParamMap &v);
	virtual ~FAMF();
	ReturnFlag run_();
	int removeOverlapping();
	double getAvgCurRadius();
private:
	void createSwarms();
	ReturnFlag increaseDiversity();
	bool checkIncreaseDiv();
	void updateIndis(unsigned );
	bool ifTerminating();
	int getNumHiber();
	void removeRedundentHiber();
protected:
	void increaseDimension();
	void decreaseDimension();
	void updateMemory();
	void wirteFile();
};

template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
FAMF<TypeIndi,TypeMainPop,TypeSubPop>::FAMF(ParamMap &v):TypeMainPop((v[param_popSize]),(v[param_evalCountFlag])), \
MultiPopulationCont<TypeSubPop>((v[param_subPopSize]),(0/*v[param_overlapDgre]*/))\
,mv_convered(),mv_indis(100),mp_single(new TypeMainPop()),\
m_preIndis((v[param_popSize])),m_nextIndis(0),m_prePops(0),\
ms_maxPops(0),m_clst(),mc_maxIndis(10000),\
mc_minIndis(10), mc_stepIndis((v[param_stepIndi])), mc_offPeak(v[param_peakOffset]), \
m_hibernatingRadius(0),m_minNumIndis(v[param_minNumPopSize]),m_optFound(GET_NUM_DIM,1,false,false,Program_Algorithm){

	for(int i=0;i<100;i++){
		mv_indis[i].push_back(0);
		mv_indis[i].push_back(0);
	}
	this->setCongergFactor((v[param_convFactor]));
	this->setConvergThreshold(1.e-9);
	this->m_name=(string)v[param_algName];
	m_clst.setMinGroupSize(m_minNumIndis);
	if (Global::msp_global->mp_problem->isProTag(DOP)) FAMFDerating::ms_enableDerating = false;
	else FAMFDerating::ms_enableDerating=true;
	if(FAMFDerating::msp_opt.get()==nullptr){
		FAMFDerating::msp_opt.reset(&m_optFound);
	}
	this->m_algPar<<"Global population size: "<<this->m_popsize<<"; subsize: "<<this->m_maxSubSize<<"; Overlap degree:"<<this->m_overlapDegree<<"; step: "<<mc_stepIndis;
}

template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
ReturnFlag FAMF<TypeIndi,TypeMainPop,TypeSubPop>::run_(){
	
	#ifdef OFEC_CONSOLE
		if(Global::msp_global->m_runId==0){
			mSingleObj::getSingleObj()->setProgrOutputFlag(true);
			if(mMultiModal::getPopInfor())
			mMultiModal::getPopInfor()->setOutProgFlag(true);
		}
	#endif // OFEC_CONSOLE

	createSwarms();
	int convergedIndis=0;
	bool firstAdjustment=true;
	ReturnFlag r_flag=Return_Normal;
	
	static	vector<double> waction(MAX_NUM_RUN,0),taction(MAX_NUM_RUN,0);
	static double lessPreIndi=0;
	static vector<vector<pair<double,double>>> action(MAX_NUM_RUN,vector<pair<double,double>>());
	while(!ifTerminating()){
		/*if(Global::ms_curProId==Global::msm_pro["DYN_CONT_MovingPeak"])
			cout<<Global::msp_global->mp_problem->getEvaluations()<<" "<<this->m_subPop.size()<<" "<<CAST_PROBLEM_DYN_CONT->getPeaksFound()<<endl;
		else	cout<<Global::msp_global->mp_problem->getEvaluations()<<" "<<this->m_subPop.size()<<" "<<CAST_PROBLEM_CONT->getGOpt().getNumGOptFound()<<" "<<m_optFound.getNumGOptFound()<<" "<<getNumHiber()<<endl;
		*/
		for(unsigned k=0;k<this->m_subPop.size();k++){		
			if(this->m_subPop[k]->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) continue;	
			int tres=1;
			if(CAST_PROBLEM_CONT->getGOpt().getNumOpt()==1&& k==this->findBestPop(0,FAMFPop<TypeIndi,TypeMainPop>::f_hiber,false)) tres=Global::g_arg[param_resource4BestPop];
			for(int ii=0;ii<tres;++ii){	
				if(this->m_subPop[k]->m_popsize>1){
					this->m_subPop[k]->updateCurRadius(true);
					double radius=this->m_subPop[k]->m_curRadius;
					r_flag=this->m_subPop[k]->evolve();
				
					this->handleReturnFlag(r_flag);
					this->handleReturnFlagAll(r_flag);
					HANDLE_RETURN_FLAG(r_flag)			

					if(this->m_subPop[k]->m_curRadius/(radius+0.00001)<0.9) this->m_subPop[k]->m_countRadiusDecreace=0;
					else this->m_subPop[k]->m_countRadiusDecreace++;

					vector<TypeIndi> par;
					for(auto &best:this->m_subPop[k]->m_best){
						par.push_back(move(TypeIndi()));
						par.back().self()=*best;
						if(this->m_subPop[k]->isStagnant()){
							r_flag=par.back().cauchyMove();
						}
						else r_flag=par.back().brwonianMove(this->m_subPop[k]->m_curRadius);
						if(r_flag!=Return_Normal) break;
					}
			
					for(auto& p:par){
						this->m_subPop[k]->updateBestArchive(p.self());			
					}
					this->handleReturnFlag(r_flag);
					this->handleReturnFlagAll(r_flag);
					HANDLE_RETURN_FLAG(r_flag)
			
					bool flag=false;				
				
					flag=this->m_subPop[k]->updateHiberFlag();
					if(this->m_subPop[k]->isHibernating()) this->m_subPop[k]->reduceRadius(CAST_PROBLEM_CONT->getDisAccuracy());
					if(flag){
						flag=m_optFound.isFound(this->m_subPop[k]->m_best[0]->data(),CAST_PROBLEM_CONT->getDisAccuracy(),Global::msp_global->mp_problem->getAccuracy());
						if(flag){
							this->measureMultiPop();
							int last=m_optFound.getNumGOptFound()-1;
							//cout<<last+1<<endl;
							m_optFound[last].setDeratingObj(last);
							for(auto &s:this->m_subPop){
								if(CAST_PROBLEM_CONT->getDistance(s->m_best[0]->data(),m_optFound[last].data(),DIS_EUCLIDEAN)<=m_optFound[last].getRadius() ){
									s->deratingFitness(mp_single);
								}
							}
							#ifdef OFEC_DEMON
							msp_buffer->updateFitnessLandsacpe_();
							#endif
						}
					}	

					#ifdef OFEC_DEMON
						vector<Algorithm*> vp;
						for(auto &it:this->m_subPop){
							vp.push_back(it.get());
						}
						msp_buffer->updateBuffer_(&vp);
					#endif

					if(r_flag==Return_Terminate) 	break;
				}
				if(this->m_subPop[k]->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) break;
			}
			if(r_flag==Return_Terminate) 	break;
         }
		if(r_flag==Return_Terminate) 	break;
		while(-1!=removeOverlapping());
		
		removeRedundentHiber();	

		if(this->getNumPops()>1){
            for(unsigned k=0;k<this->getNumPops();k++){
				if(this->m_subPop[k]->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) continue;	
				TypeSubPop *s= this->m_subPop[k].get();
				if(s->isConverged(this->m_convergThreshold)){
					//cout<<"converge"<<endl;
					s->setHiberFlag(true);
					s->reduceRadius(CAST_PROBLEM_CONT->getDisAccuracy());
					convergedIndis+=s->m_popsize;					
					// repair individuals in mv_congered that belong to the last environment for dimensional change
					if(mv_convered.size()>0){
						for(unsigned int i=0;i<mv_convered.size();i++){
							while(mv_convered[i]->getNumDim()>GET_NUM_DIM)	mv_convered[i]->decreaseDimension();
							while(mv_convered[i]->getNumDim()<GET_NUM_DIM) mv_convered[i]->increaseDimension();
						}
					}
					//cout<<"A population has converged, Hibernating"<<s->m_flag[0]<<" "<<fabs(s->m_best[0]->obj(0)-CAST_PROBLEM_CONT->getGOpt().obj(0))<<endl;
					bool removeflag=false;
					for(auto &i:s->m_best){						
						Solution<CodeVReal> gopt;
						if(CAST_PROBLEM_CONT->getGOpt().flagGloObj()) gopt=CAST_PROBLEM_CONT->getGOpt()[0];
						else gopt=(*i);
						double good=gopt.getObjDistance_((*i).data().m_obj)/gopt.getObjDistance_(Solution<CodeVReal>::getWorstSolutionSoFar().obj());
						if(CAST_PROBLEM_CONT->getGOpt().flagLoc()&&CAST_PROBLEM_CONT->getGOpt().getNumOpt()>1||good<0.1){
							bool flag=m_optFound.isFound(i->data(),CAST_PROBLEM_CONT->getDisAccuracy(),Global::msp_global->mp_problem->getAccuracy());
							if(flag){
								//m_optFound.appendOptima(ContOptimum(*i),true);
								int last=m_optFound.getNumGOptFound()-1;
								m_optFound[last].setDeratingObj(last);
								for(auto &s:this->m_subPop){
									if(CAST_PROBLEM_CONT->getDistance(s->m_best[0]->data(),m_optFound[last].data(),DIS_EUCLIDEAN)<=m_optFound[last].getRadius() ){
										s->deratingFitness(mp_single);
									}
								}
								#ifdef OFEC_DEMON
								msp_buffer->updateFitnessLandsacpe_();
								#endif	
							}
						}else{
							mv_convered.push_back(move(unique_ptr<TypeIndi>( new TypeIndi(*i))));
							removeflag=true;
						}
					}
					if(removeflag){
						this->deletePopulation(k);
						k--;
						this->m_convergedPops++;
					}
								
                }
            }
		}
		
		this->measureMultiPop();

		if(checkIncreaseDiv()){//		
			if(Global::msp_global->mp_problem->isProTag(DOP)){
				//wakeup only in dynamic environments
				for(auto &s:this->m_subPop) s->wakeup();
				m_optFound.clear();
			}

			int curPops=this->m_convergedPops+this->getNumPops();
		
			if(firstAdjustment){
				m_prePops=curPops;
				m_preIndis=convergedIndis+Global::msp_global->m_totalNumIndis;
				firstAdjustment=false;
			}
			updateIndis(curPops);

			// fuzzy logic theory to determine the total number of indis for next moment
			int indisNum=static_cast<int>(Global::msp_global->mp_normalAlg->NextNonStand(mv_indis[curPops][0],mv_indis[curPops][1]));
						
			//cout<<Global::msp_global->mp_problem->getEvaluations()<<" "<<curPops<<" "<<m_prePops<<" "<<indisNum<<" "<<Global::msp_global->m_totalNumIndis<<endl;
				
			int dif=curPops-m_prePops;
			double ratio=fabs(dif*1./mc_offPeak);
			if(ratio>=1){
				m_nextIndis=indisNum+mc_stepIndis*dif;
			}else if(ratio==0){
				m_nextIndis=indisNum;
			}else{
				double p=Global::msp_global->mp_uniformAlg->Next();
				if(p<=ratio) m_nextIndis=indisNum+mc_stepIndis*dif;
				else m_nextIndis=indisNum;
			}
			#ifdef OFEC_DEBUG_
			//for test purpurs result
			g_mutex.lock();
			taction[Global::msp_global->m_runId]+=1;
			action[Global::msp_global->m_runId].push_back(move(make_pair((double)Global::msp_global->mp_problem->getEvaluations(),(double)taction[Global::msp_global->m_runId])));
			if(Global::ms_curProId==Global::msm_pro["DYN_CONT_MovingPeak"]){ //
				if(curPops<CAST_PROBLEM_DYN_CONT->getNumberofPeak()&&m_nextIndis<indisNum||curPops>CAST_PROBLEM_DYN_CONT->getNumberofPeak()&&m_nextIndis>indisNum) waction[Global::msp_global->m_runId]+=1;
			}else{ 
				if(curPops<CAST_PROBLEM_CONT->getGOpt().getNumOpt()&&m_nextIndis<indisNum||curPops>CAST_PROBLEM_CONT->getGOpt().getNumOpt()&&m_nextIndis>indisNum) 
					waction[Global::msp_global->m_runId]+=1;
			}
			g_mutex.unlock();
			//end
			#endif
			m_prePops=curPops;

			if(m_nextIndis<=Global::msp_global->m_totalNumIndis){
				m_nextIndis=Global::msp_global->m_totalNumIndis+mc_stepIndis*2;
				#ifdef OFEC_DEBUG_
				g_mutex.lock();
				lessPreIndi+=1;
				g_mutex.unlock();
				#endif
			}
			r_flag=increaseDiversity();
			this->m_convergedPops=0;
			m_preIndis=m_nextIndis;
			if(mv_convered.size()>0) mv_convered.clear();
			convergedIndis=0;
			if(r_flag==Return_Terminate) 	break;
		}
	}
	#ifdef OFEC_DEBUG_
	#ifdef OFEC_CONSOLE
		//for debug or middle results
		wirteFile();	
		g_mutex.lock();
		string ss=Global::g_arg[param_workingDir];
		ss+="Result/";
		ss+=mSingleObj::getSingleObj()->m_fileName.str();
		ss+="WAR.txt";
		ofstream out(ss.c_str());
		double rate=0;
		double avgPS=0;
		int maxsize=0;
		int RUN=MAX_NUM_RUN;
		for(int i=0;RUN>i;i++){
			rate+=waction[Global::msp_global->m_runId]/taction[Global::msp_global->m_runId];
			avgPS+=std::accumulate(avgIndi[Global::msp_global->m_runId].begin(),avgIndi[Global::msp_global->m_runId].end(),0.)/avgIndi[Global::msp_global->m_runId].size();
			if(maxsize<action[i].size()) maxsize=action[i].size();
		}
		rate/=(int)RUN;
		avgPS/=(int)RUN;
		out<<"# "<<rate<<" "<<waction[Global::msp_global->m_runId]<<" "<<taction[Global::msp_global->m_runId]<<" "<<avgPS<<" "<<lessPreIndi/(int)MAX_NUM_RUN<<endl;
		for(int i=0;i<maxsize;++i){
			double mfes(0),mact(0);
			for(int j=0;RUN>j;j++){
				if(i<action[j].size()){
					mfes+=action[j][i].first;
					mact+=action[j][i].second;
				}else{
					mfes+=action[j].back().first;
					mact+=action[j].back().second;
				}
			}
			out<<mfes/RUN<<" "<<mact/RUN<<endl;
		}
		out.close();
		g_mutex.unlock();
		//end
	#endif
	#endif

	return Return_Terminate;

}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
 void FAMF<TypeIndi,TypeMainPop,TypeSubPop>:: createSwarms( ){
	if(this->m_popsize<m_minNumIndis) return;
	
	for(auto &i:this->m_best){
		FAMFDerating::derateFitness(*i);
	}
	for(auto &i:this->m_pop) FAMFDerating::derateFitness(i->representative());

	FAMFDerating::derateFitness(this->m_center);

	m_clst.setNormalizationFlag(false);
	m_clst.setSpace(0);
	m_clst.initialize(this->m_pop,this->m_popsize);
	m_clst.adaptiveClustering();

	int newswarms=0;

	for(int k=0;k< m_clst.getSize();k++){
		//cout<<"pop "<<k<<" size "<<m_clst.m_group[k].m_number<<endl;
		if(m_clst[k].getSize()<m_minNumIndis) { // keep populations with single indi for later
			mp_single->add(m_clst[k]);
			continue;
		}		
		TypeSubPop *s=new TypeSubPop(m_clst[k]);
		s->setConvergThreshold(this->m_convergThreshold);
		s->checkOverCrowd(this->m_maxSubSize);
		#ifdef OFEC_DEBUG_
		g_mutex.lock();
		avgIndi[Global::msp_global->m_runId].push_back(s->m_popsize);
		g_mutex.unlock();
		#endif
		this->addPopulation(*s);
		newswarms++;
		s=0;
	}
	//cout<<this->m_popsize<<" "<<newswarms<<" "<<this->m_subPop.size()<<endl;
	m_clst.clear_();
	this->remove(this->m_popsize);

	this->measureMultiPop();
}
 template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
int FAMF<TypeIndi,TypeMainPop,TypeSubPop>::removeOverlapping(){
	for(unsigned i=0;i<this->getNumPops();i++){
		for(unsigned j=i+1;j<this->getNumPops();j++){
			double dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_center);
			if(dist<this->m_subPop[i]->m_initialRadius||dist<this->m_subPop[j]->m_initialRadius){
				int c1=0,c2=0;
				for(int k=0;k<this->m_subPop[j]->m_popsize;k++){
					dist=this->m_subPop[i]->m_center.getDistance(this->m_subPop[j]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c1++;
				}
				for(int k=0;k<this->m_subPop[i]->m_popsize;k++){
					dist=this->m_subPop[j]->m_center.getDistance(this->m_subPop[i]->m_pop[k]->representative());
					if(dist<this->m_subPop[i]->m_initialRadius) c2++;
				}
				if(c1>0&&c2>0){
					int idx=-1;
					if(*this->m_subPop[i]>*this->m_subPop[j]){
						this->m_subPop[i]->add(*this->m_subPop[j]);
						if(Global::ms_curProId==Global::msm_pro["DYN_CONT_MovingPeak"]){
							if(this->m_subPop[j]->m_radiusQaulity>this->m_subPop[i]->m_radiusQaulity){
								this->m_subPop[i]->m_initialRadius=this->m_subPop[j]->m_initialRadius;
								this->m_subPop[i]->m_radiusQaulity=this->m_subPop[j]->m_radiusQaulity;
							}else if(0==this->m_subPop[i]->m_radiusQaulity){
								if(this->m_subPop[i]->m_initialRadius>this->m_subPop[j]->m_initialRadius){
									this->m_subPop[i]->m_initialRadius-=(this->m_subPop[i]->m_initialRadius-this->m_subPop[j]->m_initialRadius)*this->m_subPop[j]->m_initialRadius/this->m_subPop[i]->m_initialRadius;
								}else{
									this->m_subPop[i]->m_initialRadius+=(this->m_subPop[j]->m_initialRadius-this->m_subPop[i]->m_initialRadius)*this->m_subPop[i]->m_initialRadius/this->m_subPop[j]->m_initialRadius;
								}
							}
						}
						if(this->m_subPop[i]->m_popsize>this->m_maxSubSize){
							this->m_subPop[i]->sort(true,true);
							int *id=new int[this->m_subPop[i]->m_popsize-this->m_maxSubSize];
							for(int k=this->m_maxSubSize;k<this->m_subPop[i]->m_popsize;k++)
								id[k-this->m_maxSubSize]=this->m_subPop[i]->m_pop[this->m_subPop[i]->m_orderList[k]]->m_id;
							this->m_subPop[i]->remove(this->m_subPop[i]->m_popsize-this->m_maxSubSize,id);
							delete [] id;
							id=0;
						}
						this->deletePopulation(j);
						idx=j;

					}else{
						this->m_subPop[j]->add(*this->m_subPop[i]);
						if(Global::ms_curProId==Global::msm_pro["DYN_CONT_MovingPeak"]){
							if(this->m_subPop[j]->m_radiusQaulity<this->m_subPop[i]->m_radiusQaulity){
								this->m_subPop[j]->m_initialRadius=this->m_subPop[i]->m_initialRadius;
								this->m_subPop[j]->m_radiusQaulity=this->m_subPop[i]->m_radiusQaulity;
							}else if(0==this->m_subPop[j]->m_radiusQaulity){
								if(this->m_subPop[i]->m_initialRadius>this->m_subPop[j]->m_initialRadius){
									this->m_subPop[j]->m_initialRadius+=(this->m_subPop[i]->m_initialRadius-this->m_subPop[j]->m_initialRadius)*this->m_subPop[j]->m_initialRadius/this->m_subPop[i]->m_initialRadius;
								}else{
									this->m_subPop[j]->m_initialRadius-=(this->m_subPop[j]->m_initialRadius-this->m_subPop[i]->m_initialRadius)*this->m_subPop[i]->m_initialRadius/this->m_subPop[j]->m_initialRadius;
								}
							}
						}
						if(this->m_subPop[j]->m_popsize>this->m_maxSubSize){
							this->m_subPop[j]->sort(true,true);
							int *id=new int[this->m_subPop[j]->m_popsize-this->m_maxSubSize];
							for(int k=this->m_maxSubSize;k<this->m_subPop[j]->m_popsize;k++)
								id[k-this->m_maxSubSize]=this->m_subPop[j]->m_pop[this->m_subPop[j]->m_orderList[k]]->m_id;
							this->m_subPop[j]->remove(this->m_subPop[j]->m_popsize-this->m_maxSubSize,id);
							delete [] id;
							id=0;
						}
						this->deletePopulation(i);
						idx=i;
					}
					
					return idx;
				}
			}
		}
	}
	return -1;
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
ReturnFlag FAMF<TypeIndi,TypeMainPop,TypeSubPop>::increaseDiversity(){

	ReturnFlag rf=Return_Normal;
	if(m_nextIndis>mc_maxIndis) m_nextIndis=mc_maxIndis;
	if(m_nextIndis<mc_minIndis) m_nextIndis=mc_minIndis;
	if(m_nextIndis<=Global::msp_global->m_totalNumIndis){
		return Return_Normal;
	}
	int number=m_nextIndis-Global::msp_global->m_totalNumIndis;

	double avgRadius=this->getAvgInitialRadius();
	
	if(mv_convered.size()>0){
		// repair individuals in mv_congered that belong to the last environment for dimensional change
		for(unsigned int i=0;i<mv_convered.size();i++){
			if(mv_convered[i]->getNumDim()>GET_NUM_DIM)	mv_convered[i]->decreaseDimension();
			if(mv_convered[i]->getNumDim()<GET_NUM_DIM) mv_convered[i]->increaseDimension();
		}
		for(unsigned conv=0;conv<mv_convered.size();conv++){
			if(number<2) break;
			if ((!CAST_PROBLEM_DYN || CAST_PROBLEM_DYN&&CAST_PROBLEM_DYN->predictChange(2))/*|| (!CAST_PROBLEM_DYN_ONEPEAK || CAST_PROBLEM_DYN_ONEPEAK&&CAST_PROBLEM_DYN_ONEPEAK->predictChange(2))*/)
				rf=this->add(2,false,false);
			else
				rf=this->add(2,true,false);
			if(rf==Return_Terminate){ return rf;}

			this->m_pop[this->m_popsize-1]->initialize(*mv_convered[conv],avgRadius,this->m_popsize-1,this->m_popsize,false);
			this->m_pop[this->m_popsize-2]->initialize(*mv_convered[conv],avgRadius,this->m_popsize-2,this->m_popsize-1,false);

			number-=2;
		}
	}

	if(mp_single->getPopSize()>0){
		// repair individuals in single pop that belong to the last environment for dimensional change
		for( int i=0;i<mp_single->getPopSize();i++){
			if((*mp_single)[i]->getNumDim()>GET_NUM_DIM)	this->m_pop[i]->decreaseDimension();
			if((*mp_single)[i]->getNumDim()<GET_NUM_DIM) this->m_pop[i]->increaseDimension();
		}

		for(int i=0;i<mp_single->getPopSize();i++){
			if(number<2) break;
			if ((!CAST_PROBLEM_DYN || CAST_PROBLEM_DYN&&CAST_PROBLEM_DYN->predictChange(2)) /*|| (!CAST_PROBLEM_DYN_ONEPEAK || CAST_PROBLEM_DYN_ONEPEAK&&CAST_PROBLEM_DYN_ONEPEAK->predictChange(2))*/)
				rf=this->add(2,false,false);
			else
				rf=this->add(2,true,false);
			if(rf==Return_Terminate){ return rf;}

			this->m_pop[this->m_popsize-1]->initialize((*mp_single)[i]->self(),avgRadius,this->m_popsize-1,this->m_popsize,false);
			this->m_pop[this->m_popsize-2]->initialize((*mp_single)[i]->self(),avgRadius,this->m_popsize-2,this->m_popsize-1,false);

			number-=2;
		}
	}
	
	if(number>0){
		if ((!CAST_PROBLEM_DYN || CAST_PROBLEM_DYN&&CAST_PROBLEM_DYN->predictChange(number)) /*|| (!CAST_PROBLEM_DYN_ONEPEAK || CAST_PROBLEM_DYN_ONEPEAK&&CAST_PROBLEM_DYN_ONEPEAK->predictChange(number))*/) {
			rf=this->add(number,false,false);
		}else{
          	rf=this->add(number,true,false);
			if(rf==Return_Terminate){ return rf;}
		}
	}

	if(mv_convered.size()>0){
		this->add(mv_convered);
	}

	if(mp_single->getPopSize()>0){
		this->add(*mp_single);
		mp_single->remove(mp_single->getPopSize());
	}
	if(this->m_popsize<m_minNumIndis){
		rf=this->add(m_minNumIndis-this->m_popsize,false,false);
		//cout<<this->m_popsize<<" happen "<<m_minNumIndis<<endl;
	} 
	createSwarms();
	return rf;
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
bool FAMF<TypeIndi,TypeMainPop,TypeSubPop>::checkIncreaseDiv(){
	if(this->getNumPops()==0)
		return true;

	double avgRadius=0;
	int count=0;
	for(unsigned k=0;k<this->getNumPops();k++){
		double r=this->getAvgCurRadius();
		if(!this->m_subPop[k]->isHibernating()&&!this->m_subPop[k]->isStagnant(this->m_convFactor,r)){
			avgRadius+=this->m_subPop[k]->m_curRadius;
			count++;
		}
	}
	if(count>0)	avgRadius/=count;
	else return true;
	if(avgRadius<=this->m_convFactor*CAST_PROBLEM_CONT->getDomainSize()){//
		return true;
	}
	else return false;
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::updateIndis(unsigned curPops){
	if(curPops>=mv_indis.size()){
		unsigned size=mv_indis.size();
		mv_indis.resize(curPops+1);
		for(unsigned i=size;i<mv_indis.size();i++){
			mv_indis[i].push_back(0);
			mv_indis[i].push_back(0);
		}
	}
	mv_indis[curPops].push_back(m_preIndis);
	double mean=mv_indis[curPops][0];
	mv_indis[curPops][0]=(mean*(mv_indis[curPops].size()-3)+m_preIndis)/(mv_indis[curPops].size()-2);
	mv_indis[curPops][1]=0;
	for(unsigned k=2;k<mv_indis[curPops].size();k++) {
		mv_indis[curPops][1]+=(mv_indis[curPops][k]-mv_indis[curPops][0])*(mv_indis[curPops][k]-mv_indis[curPops][0]);
	}
	mv_indis[curPops][1]=sqrt(mv_indis[curPops][1]/(mv_indis[curPops].size()-2));

}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::increaseDimension(){
	TypeMainPop::increaseDimension();
	for(auto& i:mv_convered) i->increaseDimension();
	mp_single->increaseDimension();
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::decreaseDimension(){
	TypeMainPop::decreaseDimension();
	for(auto& i:mv_convered) i->decreaseDimension();
	mp_single->decreaseDimension();
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::updateMemory(){
	TypeMainPop::updateMemory();
	for(auto& i:mv_convered) i->updateMemory();
	if(mp_single->getPopSize()>0)	mp_single->updateMemory();
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
FAMF<TypeIndi,TypeMainPop,TypeSubPop>::~FAMF(){
	mv_convered.clear();
	mp_single.reset();
	mv_indis.clear();
	FAMFDerating::msp_opt.release();
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
bool FAMF<TypeIndi,TypeMainPop,TypeSubPop>::ifTerminating(){
	
	#ifdef OFEC_DEMON
		return Algorithm::ifTerminating();
	#endif

	#ifdef OFEC_CONSOLE
		if (Global::msp_global->mp_problem->isProTag(DOP)) return Algorithm::ifTerminating();
	if(Global::msp_global->mp_problem->m_name.find("FUN_")!=string::npos){
		if(CAST_PROBLEM_CONT->getGOpt().getNumGOptFound()==CAST_PROBLEM_CONT->getGOpt().getNumOpt()||Algorithm::ifTerminating()){
			this->measureMultiPop();	
			return true;
		}
	}
	return false;
	#endif
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
int FAMF<TypeIndi,TypeMainPop,TypeSubPop>::getNumHiber(){
	int num=0;
	for(auto &i:this->m_subPop){
		if(i->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) ++num;
	}
	return num;
}
template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::removeRedundentHiber(){

	for(int i=0;i<m_optFound.getNumGOptFound();++i){
		bool flag=false;
		for(decltype(this->m_subPop.size()) j=0;j<this->m_subPop.size();++j){
			if(!this->m_subPop[j]->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) continue;
			double dis=Global::msp_global->mp_problem->getDistance(m_optFound[i].data(),this->m_subPop[j]->m_best[0]->data(),DIS_EUCLIDEAN);
			if(dis<=CAST_PROBLEM_CONT->getDisAccuracy()){
				if(!flag) flag=true;
				else{
					this->m_subPop.erase(this->m_subPop.begin()+j);
					--j;
					//cout<<"remove redundent hib pop"<<endl;
				}
			}
		}
	}
}

template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
double FAMF<TypeIndi,TypeMainPop,TypeSubPop>::getAvgCurRadius(){
	if(this->m_subPop.size()==0) return 0;
	double r=0;
	int count=0;
	for(unsigned int j=0;j<this->m_subPop.size();j++){
		if(this->m_subPop[j]->m_flag[FAMFPop<TypeIndi,TypeMainPop>::f_hiber]) continue;
		r+=this->m_subPop[j]->m_curRadius;
		count++;
	}
	if(count>0) 	r/=count;
	return r;
}

template<typename TypeIndi, typename TypeMainPop, typename TypeSubPop>
void FAMF<TypeIndi,TypeMainPop,TypeSubPop>::wirteFile(){
		
	g_mutex.lock();	
	double total=0;
	for(int i=0;i<mv_indis.size();i++) total+=mv_indis[i].size()-2;
	static vector<double> distri;
	if(mv_indis.size()>distri.size()){
		distri.resize(mv_indis.size(),0);
	}
	for(int i=0;i<mv_indis.size();i++){
		distri[i]+=(mv_indis[i].size()-2)/total;
	}		

	string ss=Global::g_arg[param_workingDir];
	ss+="Result/";
	ss+=mSingleObj::getSingleObj()->m_fileName.str();
	ss+="distr.txt";

	/*if(Global::msp_global->mp_problem->isProTag(DOP)){
		ss<<Global::g_arg[param_workingDir]<<"Result/peaks"<<CAST_PROBLEM_DYN->getInitialNumPeaks()<<"numPeaksChange"<<CAST_PROBLEM_DYN->getFlagNumPeaksChange()<<"type"<<CAST_PROBLEM_DYN->getNumPeakChangeMode()<<"_distr.txt";
	}else{
		ss<<Global::g_arg[param_workingDir]<<"Result/"<<Global::msp_global->mp_problem->m_name<<"_Dim_"<<GET_NUM_DIM<<"_distr.txt";
	}*/

	ofstream out(ss.c_str());
	for(int i=0;i<distri.size();i++) out<<i<<" "<<distri[i]/(int)MAX_NUM_RUN<<endl;
	out.close();
	
	/*if(CAST_PROBLEM_CONT->getGOpt().flagGloObj()){
		static int num=0;
		if(num<CAST_PROBLEM_CONT->getGOpt().getNumGOptFound()){
			num=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
			stringstream ss;
			ss<<Global::g_arg[param_workingDir]<<"Result/"<<Global::msp_global->mp_problem->m_name<<"_Opt_"<<GET_NUM_DIM<<"_Dim.txt";
			ofstream out(ss.str().c_str());
			CAST_PROBLEM_CONT->getGOpt().printGOpt(out);
			out.close();
		}
	}else{
		static int num=0;
		if(num<m_optFound.getNumGOptFound()){
			num=m_optFound.getNumGOptFound();
			stringstream ss;
			ss<<Global::g_arg[param_workingDir]<<"Result/"<<Global::msp_global->mp_problem->m_name<<"_Opt_"<<GET_NUM_DIM<<"_Dim.txt";
			ofstream out(ss.str().c_str());
			m_optFound.printGOpt(out);
			out.close();
		}
	}*/
	g_mutex.unlock();
}
#endif
