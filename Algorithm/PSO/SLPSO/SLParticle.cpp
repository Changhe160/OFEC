#include "SLParticle.h"

SLParticle::SLParticle():Particle(),m_itersUnimpr(0),m_updateFre(0),m_learnRatio(0)
{
    //ctor
	mv_prog.resize(ms_numOperators);
	mv_monitor.resize(ms_numOperators);
}

SLParticle::~SLParticle()
{
}

SLParticle::SLParticle( const SLParticle& other):Particle(other)
{
    //copy ctor
    mv_prog.resize(ms_numOperators);
	mv_monitor.resize(ms_numOperators);

    m_itersUnimpr=other.m_itersUnimpr;
	m_updateFre=other.m_updateFre;
	m_learnRatio=other.m_learnRatio;

	copy(other.mv_prog.begin(),other.mv_prog.end(),mv_prog.begin());
	copy(other.mv_monitor.begin(),other.mv_monitor.end(),mv_monitor.begin());
}

SLParticle& SLParticle::operator=( const SLParticle& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator

    Particle::operator=(rhs);
	m_itersUnimpr=rhs.m_itersUnimpr;
	m_updateFre=rhs.m_updateFre;
	m_learnRatio=rhs.m_learnRatio;

    copy(rhs.mv_prog.begin(),rhs.mv_prog.end(),mv_prog.begin());
	copy(rhs.mv_monitor.begin(),rhs.mv_monitor.end(),mv_monitor.begin());

    return *this;
}

void SLParticle::setSelRatio(){
	for(int i=0;i<ms_numOperators;i++)	{
		if(m_flag)
			mv_prog[i].m_ratio=1./ms_numOperators;
		else{
			if(i==ms_numOperators-1) mv_prog[i].m_ratio=0.0;
			else mv_prog[i].m_ratio=1./(ms_numOperators-1);
		}

		mv_prog[i].m_minRatio=0.001;
		mv_prog[i].m_numSelected=0;
		mv_prog[i].m_numSuccess=0;
		mv_prog[i].m_rewards=0;
		mv_monitor[i]=mv_prog[i];
	}

}

void SLParticle::nonLearnToLearn(){
    for( int m=0;m<ms_numOperators;m++){
        mv_prog[m].initialize(ms_numOperators);
	}

	double sum=0;

	for(int m=0;m<ms_numOperators-1;m++) sum+=mv_monitor[m].m_ratio;

	int m;
	for( m=0;m<ms_numOperators-1;m++)	{
		mv_monitor[m].m_ratio=(double)(ms_numOperators-1.)/ms_numOperators*(mv_monitor[m].m_ratio/sum);

	}
	mv_monitor[m].initialize(ms_numOperators);

}

void SLParticle::learnToNonLearn(){
    double sum=0;
	for(int j=0;j<ms_numOperators-1;j++)
		sum+=mv_prog[j].m_ratio;
	for(int j=0;j<ms_numOperators-1;j++)
		mv_prog[j].m_ratio=mv_prog[j].m_ratio/sum;

	mv_prog[ms_numOperators-1].m_ratio=0;
	sum=0;
	for(int j=0;j<ms_numOperators-1;j++)
		sum+=mv_monitor[j].m_ratio;
	for(int j=0;j<ms_numOperators-1;j++)
		mv_monitor[j].m_ratio=mv_monitor[j].m_ratio/sum;
	mv_monitor[ms_numOperators-1].m_ratio=0;

}

void SLParticle::updateSelectionRatioMonitor(){

    if(m_flag)
		Progress::updateProgress(Global::msp_global.get(),mv_monitor,ms_numOperators);
	else
		Progress::updateProgress(Global::msp_global.get(),mv_monitor,ms_numOperators-1);

	for(int j=0;j<ms_numOperators;j++){
        mv_monitor[j].m_numSelected=0;
        mv_monitor[j].m_numSuccess=0;
        mv_monitor[j].m_rewards=0;
    }
}

void SLParticle::updateSelectionRatioProg(){
    if(m_flag)
		Progress::updateProgress(Global::msp_global.get(),mv_prog,ms_numOperators);
	else
		Progress::updateProgress(Global::msp_global.get(),mv_prog,ms_numOperators-1);

	for(int j=0;j<ms_numOperators;j++){
			mv_prog[j].m_numSelected=0;
			mv_prog[j].m_numSuccess=0;
			mv_prog[j].m_rewards=0;
	}

}

int SLParticle::selectOperator(){
    if(m_flag)	return Progress::getAction(Global::msp_global.get(),ms_numOperators,mv_prog);
	else	return Progress::getAction(Global::msp_global.get(),ms_numOperators-1,mv_prog);

}

ReturnFlag SLParticle::move(const Solution<CodeVReal> &lbest,double w, double c1){
	double l,u;
	for(int j=0;j<Global::msp_global->mp_problem->getNumDim();j++){
        CAST_PROBLEM_CONT->getSearchRange(l,u,j);
		m_vel[j]=w*m_vel[j]+c1*Global::msp_global->mp_uniformAlg->Next()*(lbest.data()[j]-data()[j]);
        if(m_vel[j]>m_vMax[j].m_max)	m_vel[j]=m_vMax[j].m_max;
        else if(m_vel[j]<m_vMax[j].m_min)		m_vel[j]=m_vMax[j].m_min;
		data()[j]+= m_vel[j];

	}
	self().validate();
	return evaluate();
}
