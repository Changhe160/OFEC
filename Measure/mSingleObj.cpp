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
// Created: 21 July 2011
// Last modified:
#include "mSingleObj.h"
#include "../Global/global.h"
#include "../Problem/DOP/DynamicContinuous.h"
#include "../Problem/FunctionOpt/BenchmarkFunction.h"
#include "../Algorithm/Algorithm.h"
#include "../Problem/DOP/MovingPeak.h"
#include "../Utility/TypeVar/typeVar.h"

unique_ptr<mSingleObj> mSingleObj::msp_perf(nullptr);
mSingleObj::mSingleObj(ParamMap &v):m_gOptIdx((int)MAX_NUM_RUN,0),m_numEvals2Suc(0),m_numEvals2Converge(0),m_progrOutputFlag(true),m_convgMode(PROGR_MEAN),m_absoluteErr(true),\
	m_offlineErrOverChanges(0),m_offlineErrOverChangesVar(0),m_offlineErrOverRuns(0), m_offlineErrOverRunsVar(0),m_meanOverRuns(0),m_meanOverRunsVar(0),\
	m_bestSoFar((int)MAX_NUM_RUN,0),m_avgEvals(0),m_avgCevals(0),m_avgTevals(0){
    //ctor
	int maxNumrun=MAX_NUM_RUN;
	if (Global::msp_global->mp_problem->isProTag(DOP))
		m_changeFre=int(v[param_changeFre]);
	mpp_data.resize(maxNumrun);
	mpp_convergeTime.resize(maxNumrun); 
	m_duration.resize(maxNumrun);

	if((v[param_gOptFlag])==true){
		mpp_gOpt.resize(maxNumrun);
    }
}

mSingleObj::~mSingleObj(){
 
}

mSingleObj& mSingleObj::operator=(const mSingleObj& rhs){
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

double mSingleObj::getConvergeSpeed2Best(){
    return m_speed2Best;
}
double mSingleObj::getConvergeSpeed2GOpt(){
    return m_speed2gOpt;

}
double mSingleObj::getConvergeSpeed(){

    return m_speed;

}
void mSingleObj::calculateConvergenceTime(){
//track back from the end of mpp_data[i], the convergence time of mpp_convergeTime[i] is the index of mpp_data[i] when the value changes.
//e.g., mpp_data has ten elements: 10 9 8 8 7 6 3 2 2 2. the convergence time is 7 as  mpp_data[7] is different from it's previous one mpp_data[6]
    m_numEvals2Converge=0;
	int maxNumrun=MAX_NUM_RUN;
    for(int i=0;maxNumrun>i;i++) {
        for(int j=m_numChanges-1;j>=0;j--){
            for(int k=m_numRecords-(m_numChanges-j-1)*m_recordsPerChange-1;k>=m_numRecords-(m_numChanges-j)*m_recordsPerChange;k--){
                if(k==m_numRecords-(m_numChanges-j)*m_recordsPerChange){
                    mpp_convergeTime[i][j]=k;
                    m_numEvals2Converge+=(k-(m_numRecords-(m_numChanges-j)*m_recordsPerChange)+1)*Global::g_arg[param_sampleFre];
                    break;
                }
                if(mpp_data[i][k]!=mpp_data[i][k-1]) {
                    mpp_convergeTime[i][j]=k;
                    m_numEvals2Converge+=(k-(m_numRecords-(m_numChanges-j)*m_recordsPerChange)+1)*Global::g_arg[param_sampleFre];
                    break;
                }
            }
        }
    }
    //the average number of evaluations before the algorithm converges
    m_numEvals2Converge/=(maxNumrun*m_numChanges);

}
void mSingleObj::getStatInfor(double &mean, double &var, double &worst, double &best){
    mean=m_mean; var=m_var;
	if(mpp_gOpt.size()>0){
		worst=m_largest;
        best=m_smallest;
	}else{
		if(m_comp==MIN_OPT){
			worst=m_largest;
			best=m_smallest;
		}else{
			best=m_largest;
			worst=m_smallest;

		}
	}
}
double mSingleObj::getSucRate(){
    return m_sucRate;

}
double mSingleObj::getPerformance(){
    return m_performance;
}
mSingleObj* mSingleObj::getSingleObj(){
	 return mSingleObj::msp_perf.get();
}
void mSingleObj::initialize(ParamMap &v){

    if(mSingleObj::msp_perf)        return;
    mSingleObj::msp_perf=unique_ptr<mSingleObj>(new mSingleObj(v));
}

void mSingleObj::record(Global *glob,double rObj){
	int evals=glob->mp_problem->getEvaluations();
	if(Global::msp_global->mp_algorithm!=nullptr&&Global::msp_global->mp_algorithm->ifTerminating()){
		m_avgEvals+=evals;
		m_avgCevals+=glob->mp_problem->cevals();
		m_avgTevals+=glob->mp_problem->getTevals();
	}
	if(Global::msp_global->mp_algorithm!=nullptr&&Global::msp_global->mp_algorithm->ifTerminated()) return;

	if (glob->mp_problem->isProTag(DOP))
	{
		if(evals%(m_recordsPerChange*Global::g_arg[param_sampleFre])==1) 
			m_bestSoFar[glob->m_runId]=rObj;
	}
	else
	{
		if(evals==1)
			m_bestSoFar[glob->m_runId]=rObj;
	}

	if(glob->mp_problem->getOptType()==MIN_OPT){
        if(m_bestSoFar[glob->m_runId]>rObj) m_bestSoFar[glob->m_runId]=rObj;
    }else{
       if(m_bestSoFar[glob->m_runId]<rObj) m_bestSoFar[glob->m_runId]=rObj;
    }

	if(evals%Global::g_arg[param_sampleFre]==0||glob->mp_algorithm!=nullptr&&glob->mp_algorithm->ifTerminating()){
		mpp_data[glob->m_runId].push_back(m_bestSoFar[glob->m_runId]);
        //cout<<glob->mp_problem->getEvaluations()/Global::g_arg[param_sampleFre]-1<<"     "
        //<<mpp_data[Global::g_runIdx][glob->mp_problem->getEvaluations()/Global::g_arg[param_sampleFre]-1]<<endl;
    }

}
void  mSingleObj::calculatePerformance(){

    calculateConvergenceTime();
	
	int maxNumrun=MAX_NUM_RUN;
    //fitness decrease/increase per evaluation: m_speed=f(0)-f(covergeTime)/covergeTime
    m_speed=0;
    for(int i=0;maxNumrun>i;i++){
         for(int j=m_numChanges-1;j>=0;j--) {
			 int k=(mpp_convergeTime[i][j]/m_recordsPerChange)*m_recordsPerChange;
                 m_speed+=fabs(mpp_data[i][k]-mpp_data[i][mpp_convergeTime[i][j]])
                 /(Global::g_arg[param_sampleFre]*(1+mpp_convergeTime[i][j]-k));
         }
    }
    m_speed/=(maxNumrun*m_numChanges);
    //the convergence speed to the best result achieved by the algorithm
    m_speed2Best=0;
    for(int i=0;maxNumrun>i;i++){
         for(int j=0;j<m_numChanges;j++){
            double t=0;
            for(int k=m_recordsPerChange*j;k<=mpp_convergeTime[i][j];k++){
                t+=fabs(mpp_data[i][k]-mpp_data[i][mpp_convergeTime[i][j]]);
            }
            m_speed2Best+=t/((mpp_convergeTime[i][j]-m_recordsPerChange*j+1)*Global::g_arg[param_sampleFre]);
         }

    }
    m_speed2Best/=(maxNumrun*m_numChanges);
    //the convergence speed to the global optimum
    if(mpp_gOpt.size()>0){
        m_speed2gOpt=0;
        for(int i=0;maxNumrun>i;i++){
             for(int j=0;j<m_numChanges;j++){
                double t=0;
                for(int k=m_recordsPerChange*j;k<=mpp_convergeTime[i][j];k++){
                    t+=fabs(mpp_data[i][k]-mpp_gOpt[i][j]);
                }
                m_speed2gOpt+=t/((mpp_convergeTime[i][j]-m_recordsPerChange*j+1)*Global::g_arg[param_sampleFre]);
             }

        }
        m_speed2gOpt/=(maxNumrun*m_numChanges);
    }

    // the mean, smallest and largest errors to the global optimum
    if(mpp_gOpt.size()>0){
        m_largest=m_smallest=fabs(mpp_data[0][m_recordsPerChange-1]-mpp_gOpt[0][0]);
        m_mean=0;
        vector<double> mean(m_numChanges);

        for(int j=1;j<=m_numChanges;j++){
            double er;
             mean[j-1]=0;
            for(int i=0;maxNumrun>i;i++){
                er=fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1]);
                if(m_largest<er) {
					m_largest=er;
				}
                if(m_smallest>er) m_smallest=er;
                mean[j-1]+=er;
            }
            mean[j-1]/=maxNumrun;
            m_mean+=mean[j-1];
        }
        m_mean/=m_numChanges;
        m_var=0;

        for(int j=1;j<=m_numChanges;j++){
             long double var=0;
            for(int i=0;maxNumrun>i;i++){
                long double x=fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1])-mean[j-1];
                 var+=x*x;
            }

            var=sqrt(var/(maxNumrun));
            m_var+=var;
        }
        m_var/=m_numChanges;
      

		//****03.08.2013***** average over runs *************//
		m_meanOverRuns=0;
		m_meanOverRunsVar=0;

		mean.resize(maxNumrun);
		for(int i=0;maxNumrun>i;i++){
			mean[i]=0;
			for(int j=1;j<=m_numChanges;j++){
				mean[i]+=fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1]);
			}
			mean[i]/=m_numChanges;
			m_meanOverRuns+=mean[i];
		}
		m_meanOverRuns/=maxNumrun;
		for(int i=0;maxNumrun>i;i++){
			m_meanOverRunsVar+=(mean[i]-m_meanOverRuns)*(mean[i]-m_meanOverRuns);
		}
		m_meanOverRunsVar=sqrt(m_meanOverRunsVar/maxNumrun);		
    }else{
        m_largest=m_smallest=mpp_data[0][m_recordsPerChange-1];
        m_mean=0;
        double *mean=new double[m_numChanges];

        for(int j=1;j<=m_numChanges;j++){
            double er;
             mean[j-1]=0;
            for(int i=0;maxNumrun>i;i++){
                er=mpp_data[i][j*m_recordsPerChange-1];
                if(m_largest<er) m_largest=er;
                if(m_smallest>er) m_smallest=er;
                mean[j-1]+=er;
            }
            mean[j-1]/=maxNumrun;
            m_mean+=mean[j-1];
        }
        m_mean/=m_numChanges;
        m_var=0;

        for(int j=1;j<=m_numChanges;j++){
             long double var=0;
            for(int i=0;maxNumrun>i;i++){
                 long double x=mpp_data[i][j*m_recordsPerChange-1] -mean[j-1];
                 var=var+x*x;
            }

            var=sqrt(var/(maxNumrun));
            m_var+=var;
        }
        m_var/=m_numChanges;
        delete [] mean;

    }

    //the success rate and the average number of evaluations to achieve the given accuracy level
    if(mpp_gOpt.size()>0){
        m_sucRate=0;
        int numSuc=0;
        m_numEvals2Suc=0;
		if(!gIsDynamicProlem()){
           for(int i=0;maxNumrun>i;i++){
                 for(int j=1;j<=m_numChanges;j++){
                    if(fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1])<= m_accuracy) {
                        m_sucRate+=1;
                        numSuc++;
                        for(int k=m_recordsPerChange*(j-1);k<=mpp_convergeTime[i][j-1];k++){
                            if(fabs(mpp_data[i][k]-mpp_gOpt[i][j-1])<m_accuracy){
                                m_numEvals2Suc+=(k-m_recordsPerChange*(j-1)+1);
                                break;
                            }
                        }
                    }

                 }
            }
		}else{
            for(int i=0;maxNumrun>i;i++){
                 for(int j=1;j<=m_numChanges;j++){                 
                     if(fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1])<= m_accuracy) {
                        m_sucRate+=1;
                        numSuc++;
                        for(int k=m_recordsPerChange*(j-1);k<=mpp_convergeTime[i][j-1];k++){
                            if(fabs(mpp_data[i][k]-mpp_gOpt[i][j-1])<m_accuracy){
                                m_numEvals2Suc+=(k-m_recordsPerChange*(j-1)+1);
                                break;
                            }
                        }
                    }
                 }
            }

        }
        m_sucRate/=(maxNumrun*m_numChanges);
        if(numSuc>0)        m_numEvals2Suc=m_numEvals2Suc*Global::g_arg[param_sampleFre]/numSuc;

    }

    // performance of the algorithm according to its convergene speed to GOPT and the best results obtained
     if(mpp_gOpt.size()>0){
        m_performance=0;
		m_perfVar=0;
		double *mean_per=new double[m_numChanges];
		double *perf=new double[maxNumrun];
		stringstream ss;
		ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Perf.txt";
		ofstream out(ss.str().c_str());
		//NOTE: calculate average score: for each change inerval, average the score over all runs, and then calculate the 
		//average score over all changes.  

		for(int j=0;j<m_numChanges;j++){
			mean_per[j]=0;
			for(int i=0;maxNumrun>i;i++){
                double t=0,gap;
                if(m_comp==MIN_OPT)  gap=fabs(mpp_data[i][m_recordsPerChange*(j+1)-1])+1;
                else gap=fabs(mpp_data[i][m_recordsPerChange*j])+1;

                for(int k=m_recordsPerChange*j;k<m_recordsPerChange*(j+1);k++){
                    if(m_comp==MIN_OPT){
                        ///Note: preprocess the data to make sure they are greater than 0 to avoid overflow error
                        t+=1-(mpp_gOpt[i][j]+gap)/(mpp_data[i][k]+gap);
                    }else{
                        t+=1-(mpp_data[i][k]+gap)/(mpp_gOpt[i][j]+gap);
                    }

                }
                 t/=m_recordsPerChange;
                double best;
                int bst_Idx=m_recordsPerChange*(j+1)-1;

                if(m_comp==MIN_OPT){
                    best=(mpp_gOpt[i][j]+gap)/(mpp_data[i][bst_Idx]+gap);
                }else{
                    best=(mpp_data[i][bst_Idx]+gap)/(mpp_gOpt[i][j]+gap);
                }
				best=best/(1+t);
                m_performance+=best;
				mean_per[j]+=best;
				perf[i]=best;
             }
			double var=0;
			 mean_per[j]/=maxNumrun;
			 for(int i=0;maxNumrun>i;i++) var+=(perf[i]-mean_per[j])*(perf[i]-mean_per[j]);
			 var=sqrt(var/maxNumrun);
			 m_perfVar+=var;
		
			 out<<mean_per[j]<<" ";

        }
		out.close();
        m_performance/=(maxNumrun*m_numChanges);
		m_perfVar/=m_numChanges;
		delete [] mean_per;
		delete [] perf;
    }

	// calucate offline error over changes
	 if(mpp_gOpt.size()>0){
		 m_offlineErrOverChanges=0;
		 m_offlineErrOverChangesVar=0;
		double *mean=new double[m_numChanges];
		double *perf=new double[maxNumrun];
		
		for(int j=0;j<m_numChanges;j++){
			mean[j]=0;
			for(int i=0;i<maxNumrun;i++){
                double t=0;
                for(int k=m_recordsPerChange*j;k<m_recordsPerChange*(j+1);k++){
                    t+=fabs(mpp_gOpt[i][j]-mpp_data[i][k]);
                }
                 t/=m_recordsPerChange;
                
                m_offlineErrOverChanges+=t;
				mean[j]+=t;
				perf[i]=t;
             }
			
			 mean[j]/=maxNumrun;
			 double var=0;
			 for(int i=0;maxNumrun>i;i++) var+=(perf[i]-mean[j])*(perf[i]-mean[j]);
			 var=sqrt(var/maxNumrun);
			 m_offlineErrOverChangesVar+=var;
        }
		
        m_offlineErrOverChanges/=(maxNumrun*m_numChanges);
		m_offlineErrOverChangesVar/=m_numChanges;
		
		delete [] mean;
		delete [] perf;
    }

	 // calucate offline error over runs with traditional ways
	 if(mpp_gOpt.size()>0){
		 m_offlineErrOverRuns=0;
		 m_offlineErrOverRunsVar=0;
		
		double *perf=new double[maxNumrun];
		
		for(int i=0;maxNumrun>i;i++){
			perf[i]=0;
			for(int j=0;j<m_numChanges;j++){ 
                for(int k=m_recordsPerChange*j;k<m_recordsPerChange*(j+1);k++){
                    perf[i]+=fabs(mpp_gOpt[i][j]-mpp_data[i][k]);
                }
             }
			perf[i]/=(m_recordsPerChange*m_numChanges);
	 
        }
		
		for(int i=0;maxNumrun>i;i++){
			m_offlineErrOverRuns+=perf[i];
		}
		m_offlineErrOverRuns/=maxNumrun;
		for(int i=0;maxNumrun>i;i++){
			m_offlineErrOverRunsVar+=(perf[i]-m_offlineErrOverRuns)*(perf[i]-m_offlineErrOverRuns);	 
		}
		m_offlineErrOverRunsVar=sqrt(m_offlineErrOverRunsVar/maxNumrun);

		delete [] perf;
    }

	//2016.6.21
	 
	 for (auto &i : m_duration) {
		 m_meanDuration += i.count();
	 }
	 m_meanDuration /= m_duration.size();
}
void mSingleObj::addGOpt(int runId,double rGOptObj){
	mpp_gOpt[runId].push_back(rGOptObj);
}
void mSingleObj::resetGOptIndx(){
    for(auto& i:m_gOptIdx) i=0;
}
void mSingleObj::setFileName(stringstream &rName){
    m_fileName.str(rName.str());
}
void mSingleObj::setFileName(ParamMap &v){
	m_fileName.str("");
	for(auto &i:v){
		for(auto &j:Global::msm_param){
			if(i.first==param_gOptFlag||i.first==param_algId||i.first==param_proId||i.first==param_flagNoise||\
				i.first==param_flagNumPeakChange||i.first==param_flagTimeLinkage||i.first==param_numRun||\
				i.first==param_numTask||i.first==param_minNumPopSize||i.first==param_hibernatingRadius||\
				i.first==param_solutionValidationMode||i.first==param_evalCountFlag||\
				i.first==param_workingDir||i.first==param_sampleFre||i.first==param_maxEvals||i.first==param_flagNumPeakChange||\
				i.first==param_peakNumChangeMode||i.first==param_dataDirectory1) continue;
			if(i.first==j.second){			
				m_fileName<<j.first.substr(6)<<i.second<<"_";			
				break;
			}
		}		
	}
}
void mSingleObj::setProgrOutputFlag(bool rFlag){
    m_progrOutputFlag=rFlag;
}
void mSingleObj::setConvgProgMode(ConvgProgeMode rMode){
   m_convgMode=rMode;
}
void mSingleObj::setAlgParameter(stringstream & rPar){
    m_algPar.str(rPar.str());
}
void mSingleObj::outputResult(){
	calculateKeyParam();
    calculatePerformance();
	int maxNumrun=MAX_NUM_RUN;
    if(0){//m_progrOutputFlag
        stringstream ss;
        ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Progr.txt";
        ofstream out(ss.str().c_str());
        double *data=new double[maxNumrun];
		vector<int> idx;
		stringstream oss;
        for(int i=0;i<m_numRecords;){
            for(int j=0;maxNumrun>j;j++){
				if(m_absoluteErr){
					// absolute error(offline error)
					if(mpp_gOpt.size()>0) data[j]=fabs(mpp_data[j][i]-mpp_gOpt[j][i/m_recordsPerChange]);
					else 	data[j]=mpp_data[j][i];
				}else{
					double gap;
					if(m_comp==MIN_OPT)  gap=fabs(mpp_data[j][m_recordsPerChange*(i+1)-1])+1;
					else gap=fabs(mpp_data[j][m_recordsPerChange*i])+1;

					if(m_comp==MIN_OPT){
						///Note: preprocess the data to make sure they are greater than 0 to avoid overflow error
						data[j]=1-(mpp_gOpt[j][i/m_recordsPerChange]+gap)/(mpp_data[j][i]+gap);
					}else{
						data[j]=1-(mpp_data[j][i]+gap)/(mpp_gOpt[j][i/m_recordsPerChange]+gap);
					}
				}

            }
          
			gQuickSort<double*>(data,maxNumrun,idx);
            double x;
            if(m_convgMode==PROGR_MEAN){
                x=0;
                for(int j=0;maxNumrun>j;j++) x+=data[idx[j]];
                oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<x/maxNumrun<<endl;
            }else if(m_convgMode==PROGR_MEDEAN){
                 oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<data[idx[maxNumrun/2]]<<endl;
            }else if(m_convgMode==PROGR_WORST){
                if(m_comp==MIN_OPT) oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<data[idx[maxNumrun-1]]<<endl;
                else  oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<data[idx[0]]<<endl;
            }else if(m_convgMode==PROGR_BEST){
                if(m_comp==MIN_OPT) oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<data[idx[0]]<<endl;
				else  oss<<(i+1)*Global::g_arg[param_sampleFre]<<" " <<data[idx[maxNumrun-1]]<<endl;
            }else{
                //...
            }
			i++;
			try{
				if(i%2000==0){
					out<<oss.str();
					oss.str("");
				}
			}
			catch(...){
				cout<<"memory allocation failer"<<endl;
			}
			
        }
        delete []data;
        out.close();
    }

    stringstream ss;
	ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Sta.txt";
    ofstream out(ss.str().c_str());
	out<<"Algorithm: "<<gGetAlgorithmName(Global::ms_curAlgId)<<endl;
    out<<"Algorithm Parameters: "<<m_algPar.str()<<endl;
    out<<"Number of runs: "<<maxNumrun<<endl;
    out<<"Problem: "<<gGetProblemName(Global::ms_curProId)<<endl;
    out<<"Problem Parameters: "<<m_proPar.str()<<endl;
 
    double mean, var, best, worst;
    getStatInfor(mean,var,worst,best);
	out<<"Mean of Best&STD over Runs: "<<m_meanOverRuns<<"+"<<m_meanOverRunsVar<<endl;
    out<<"Best: "<<best<<endl;
    out<<"Worst: "<<worst<<endl;
	out<<"Total evals: "<<m_avgTevals/maxNumrun<<endl;
	out<<"Countable evals: "<<m_avgCevals/maxNumrun<<endl;
	out<<"Useful evals: "<<m_avgEvals/maxNumrun<<endl;
	out << "Elapsed time(s): " << m_meanDuration << endl;
    if(mpp_gOpt.size()>0){
       //out<<"ConvergeSpeed2GOpt: "<<m_speed2gOpt<<endl;
       // out<<"Sucess rate: "<<m_sucRate<<endl;
        out<<"Performance&STD: "<<m_performance<<"+"<<m_perfVar<<endl;
		out<<"Offline error over runs&STD: "<<m_offlineErrOverRuns<<"+"<<m_offlineErrOverRunsVar<<endl;
    }
	out.close();
	
	ss.str("");
	ss<<Global::g_arg[param_workingDir]<<"Result/"<<m_fileName.str()<<"Err.txt";
	out.open(ss.str().c_str());
	
	for(int i=0;maxNumrun>i;i++){
		double err=0;
        for(int j=1;j<=m_numChanges;j++){
			if(mpp_gOpt.size()>0)
            err+=fabs(mpp_data[i][j*m_recordsPerChange-1]-mpp_gOpt[i][j-1]);
			else err+=mpp_data[i][j*m_recordsPerChange-1];
        }
		out<<err/m_numChanges<<endl;
    }
	out.close();
}

double mSingleObj::getBestSoFar(int runid){
	return m_bestSoFar[runid];

}
void mSingleObj::setAbsoluteErrFlag(bool flag){
	m_absoluteErr=flag;
}

void mSingleObj::setAccuracy(double acc){
	m_accuracy=acc;
}
void mSingleObj::setCompareType(Compare comp){
	m_comp=comp;
}
void mSingleObj::setProParameter(stringstream & rPar){
	m_proPar.str(rPar.str());
}

void mSingleObj::calculateNumRecords()
{
	int maxRun=MAX_NUM_RUN;
	int max=mpp_data[0].size();
	for(int i=1;i<maxRun;i++)
		if(max<mpp_data[i].size())
			max=mpp_data[i].size();
	for(int i=0;i<maxRun;i++)
	{
		int size=mpp_data[i].size();
		double temp=mpp_data[i][size-1];
		for(int j=size;j<max;j++)
			mpp_data[i].push_back(temp);
	}
	m_numRecords=max;
}


void mSingleObj::calculateKeyParam()
{
	calculateNumRecords();
	int maxRun=MAX_NUM_RUN;
	m_numChanges=1;
	m_recordsPerChange=m_numRecords;
	for(int i=0;maxRun>i;i++) mpp_convergeTime[i].resize(m_numChanges);
}

void mSingleObj::deleteSingleObj(){
	msp_perf.reset();
}

void mSingleObj::setDuration(const chrono::duration<double> & et, int runid) {
	m_duration[runid ]= et;
}