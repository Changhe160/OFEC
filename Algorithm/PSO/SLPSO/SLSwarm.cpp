#include "SLSwarm.h"
#ifdef OFEC_DEMON
#include "../../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif

int SLPSO::ms_updateFre=-1;
double SLPSO::ms_learnRatio=-1;
float SLPSO::ms_ratioLearnToGbest=-1;


SLPSO::SLPSO(ParamMap &v): Swarm<CodeVReal,SLParticle>(v[param_popSize],true)
{
    //ctor
	if(Global::g_arg.find(param_maxEvals)==Global::g_arg.end())
		Global::g_arg[param_maxEvals]=10000*Global::msp_global->mp_problem->getNumDim();
	initializePara(0.89f,1.496f,1.496f,true,0.89f,0.4f);
    setParameters();
}


void SLPSO::updateParameters()
{
	vector<int> rindex(m_popsize);
	Global::msp_global->initializeRandomArray<vector<int> >(rindex,m_popsize);

	if(ms_learnRatio<0)
		for(int i=0;i<m_popsize;i++){
			m_pop[rindex[i]]->m_learnRatio=(1-exp(-pow(i*1.6/m_popsize,4)))>0.05?(1-exp(-pow(i*1.6/m_popsize,4))):0.05;
		}
	else
		for(int i=0;i<m_popsize;i++){
			m_pop[i]->m_learnRatio=ms_learnRatio;
		}

	Global::msp_global->initializeRandomArray<vector<int> >(rindex,m_popsize);

	if(ms_updateFre==-1)
		for(int i=0;i<m_popsize;i++){
			m_pop[rindex[i]]->m_updateFre=10*exp(-pow(i*1.6/m_popsize,4))>1?10*exp(-pow(i*1.6/m_popsize,4)):1;
		}
	else
		for(int i=0;i<m_popsize;i++){
			m_pop[i]->m_updateFre=ms_updateFre;
		}

	updateLearnToGbest();
}

void SLPSO::setParameters(){

    for(int i=0;i<m_popsize;i++){
        if(ms_learnRatio<0)
        m_pop[i]->m_learnRatio=(1-exp(-pow(i*1.6/m_popsize,4)))>0.05?(1-exp(-pow(i*1.6/m_popsize,4))):0.05;
        else
            m_pop[i]->m_learnRatio=ms_learnRatio;

        if(ms_updateFre==-1)
        m_pop[i]->m_updateFre=10*exp(-pow(i*1.6/m_popsize,4))>1?10*exp(-pow(i*1.6/m_popsize,4)):1;
        else
            m_pop[i]->m_updateFre=ms_updateFre;
    }

    m_numLearnToGbest=m_popsize*ms_ratioLearnToGbest;
    MRandToBest(m_numLearnToGbest);
    for(int i=0;i<m_popsize;i++) m_pop[i]->setSelRatio();
}

void SLPSO::calculateNumLearning(const int sfes){
	if(m_popsize<=0) return;
	if(ms_ratioLearnToGbest<0)
		m_numLearnToGbest=(int)(m_popsize*(1-exp(-100*pow((double)(Global::msp_global->mp_problem->getEvaluations()-sfes)/(Global::g_arg[param_maxEvals]-sfes),3))));	
}


void SLPSO::updateLearnToGbest(){

	int PrenumLearnToGbest=0,CurnumLearnToGbest=0;
	for(int i=0;i<m_popsize;i++) if(m_pop[i]->getFlag()) PrenumLearnToGbest++;

	vector<bool> preLearnToGbestInfor(m_popsize);
	for(int i=0;i<m_popsize;i++) preLearnToGbestInfor[i]=m_pop[i]->getFlag();
	vector<int> preLearnToGbest;
	int k=0;
	if(PrenumLearnToGbest>0){
		preLearnToGbest.resize(PrenumLearnToGbest);
		 k=0;
		for(int i=0;i<m_popsize;i++){
			if(m_pop[i]->getFlag()) preLearnToGbest[k++]=i;
			if(k>=PrenumLearnToGbest) break;
		}
	}

	MRandToBest(m_numLearnToGbest);

	for(int i=0;i<m_popsize;i++) if(m_pop[i]->getFlag()) CurnumLearnToGbest++;
	vector<int> curLearnToGbest;

	if(CurnumLearnToGbest>0){
		curLearnToGbest.resize(CurnumLearnToGbest);
		k=0;
		for(int i=0;i<m_popsize;i++){
			if(m_pop[i]->getFlag()) curLearnToGbest[k++]=i;
			if(k>=CurnumLearnToGbest) break;
		}
	}


	for(int i=0;i<PrenumLearnToGbest;i++){
		if(!m_pop[preLearnToGbest[i]]->getFlag()) 		m_pop[preLearnToGbest[i]]->learnToNonLearn();

	}
	for(int i=0;i<CurnumLearnToGbest;i++){
		if(!preLearnToGbestInfor[curLearnToGbest[i]])	m_pop[curLearnToGbest[i]]->nonLearnToLearn();

	}

}

ReturnFlag SLPSO::evolve(){

	ReturnFlag rf=Return_Normal;
	if(m_popsize<=0)	return Return_Terminate;

	Solution<CodeVReal>  t_p;
	int numDim=Global::msp_global->mp_problem->getNumDim();
	double *avg_v=new double[numDim];

	for(int j=0;j<numDim;j++){
		avg_v[j]=0;
		for(int i=0;i<m_popsize;i++){
			avg_v[j]+=fabs(m_pop[i]->getVel()[j]);
		}
		avg_v[j]/=m_popsize;
	}
	int k;
	vector<int> rindex(m_popsize);
	Global::msp_global->initializeRandomArray<vector<int> >(rindex,m_popsize);

	for(int j=0;j<m_popsize;j++)
	{
		int i=rindex[j];
		k=i;
		if(m_pop[i]->m_itersUnimpr>m_pop[i]->m_updateFre){
			m_pop[i]->updateSelectionRatioProg();
			m_pop[i]->m_itersUnimpr=0;
		}
		int sel=m_pop[i]->selectOperator();
		t_p=m_pop[i]->self();
		int nearest;
		switch (sel+1){
		case 1:
			rf=m_pop[i]->move(m_pop[i]->representative(),m_W,m_C1);
			if(rf!=Return_Normal) 
			{
				delete [] avg_v;
				avg_v=0;
				return rf;
			}
			break;
		case 2:
			rf=m_pop[i]->NormalMutation(avg_v);
			if(rf!=Return_Normal) 
			{
				delete [] avg_v;
				avg_v=0;
				return rf;
			}
			break;
		case 3:
			nearest=Global::msp_global->getRandInt(0,m_popsize);
			if(m_popsize>1){
				while(nearest==i) nearest=Global::msp_global->getRandInt(0,m_popsize);
			}
			if(m_pop[i]->representative()>m_pop[nearest]->representative()){
				i=nearest;
				t_p=m_pop[i]->self();
				rf=m_pop[i]->move(m_pop[k]->representative(),m_W,m_C1);
				if(rf!=Return_Normal) 
				{
					delete [] avg_v;
					avg_v=0;
					return rf;
				}
			}else{
				rf=m_pop[i]->move(m_pop[nearest]->representative(),m_W,m_C1);
				if(rf!=Return_Normal) 
				{
					delete [] avg_v;
					avg_v=0;
					return rf;
				}
			}
			break;
		case 4:
			rf=m_pop[i]->move(*m_best[0],m_W,m_C1);
			if(rf!=Return_Normal) 
			{
				delete [] avg_v;
				avg_v=0;
				return rf;
			}
			break;
		default:
			delete [] avg_v;
			avg_v=0;
			throw myException("operator selection error @SLPSO::evolve()");
			break;
		}

		if(m_pop[i]->self()>(m_pop[i]->representative())){
			m_pop[i]->representative()=m_pop[i]->self();
			if(m_pop[i]->self()>(*m_best[0]))
				*m_best[0]=*m_pop[i];
		}

		m_pop[i]->mv_prog[sel].m_numSelected++;

		if(m_pop[i]->self()>(t_p)){
			m_pop[i]->mv_prog[sel].m_numSuccess++;
			m_pop[i]->mv_prog[sel].m_rewards+=fabs(t_p.obj(0)-m_pop[i]->obj(0));
		}else
			m_pop[i]->m_itersUnimpr++;

		if(m_pop[i]->self()>(t_p))
		{
			rf=updateBest(i,m_pop[i]->m_learnRatio);
			if(rf!=Return_Normal) 
			{
				delete [] avg_v;
				avg_v=0;
				return rf;
			}
		}
	}

	updateParameters();
	m_iter++;

	delete [] avg_v;
	avg_v=0;

	return rf;
}


int SLPSO::getNumLearnTogbest(){

	int n=0;
	for(int i=0;i<m_popsize;i++) if(m_pop[i]->getFlag()) n++;

	return n;
}

void SLPSO::MRandToBest( int num){

	if(num<=0){
		for(int i=0;i<m_popsize;i++)	m_pop[i]->setFlag(0);
		return;
	}
	if(num>m_popsize) num=m_popsize;

	for(int i=0;i<m_popsize;i++)	m_pop[i]->setFlag(0);
	vector<int> l;
	for(int i=0;i<m_popsize;i++) l.push_back(i);

	for(int j=0;j<num;j++){
		int k=Global::msp_global->getRandInt(0,l.size()-1);
		m_pop[l[k]]->setFlag(1);
		l.erase(l.begin()+k);
	}
	l.clear();
}


ReturnFlag SLPSO::run_()
{
	ReturnFlag rf=Return_Normal;

	 #ifdef OFEC_CONSOLE
	if(Global::msp_global->m_runId==0){
		mSingleObj::getSingleObj()->setProgrOutputFlag(true);
		if(mMultiModal::getPopInfor())
		mMultiModal::getPopInfor()->setOutProgFlag(true);
	}
	#endif // OFEC_CONSOLE

	while(!ifTerminating())
	{
        calculateNumLearning(0);
        rf=evolve();
		#ifdef OFEC_CONSOLE
		if(mMultiModal::getPopInfor()){
			int peaksf;		
			peaksf=CAST_PROBLEM_CONT->getGOpt().getNumGOptFound();
			mMultiModal::getPopInfor()->input(Global::msp_global.get(), Global::msp_global->mp_problem->getEvaluations(),\
				Global::msp_global->m_totalNumIndis,1,peaksf,\
				CAST_PROBLEM_CONT->getGOpt().getNumOpt(),0,0,0,0,\
				0,0,CAST_PROBLEM_CONT->getGOpt().isAllFound());
		}
		#endif
		
		#ifdef OFEC_DEMON
			vector<Algorithm*> vp;	
			vp.push_back(this);	
			msp_buffer->updateBuffer_(&vp);
		#endif

		if(rf==Return_Terminate) break;

		int maxEvals=Global::g_arg[param_maxEvals];
		initializePara(m_maxW-(m_maxW-m_minW)*Global::msp_global->mp_problem->getEvaluations()/maxEvals,m_C1,m_C2,false,m_maxW,m_minW);
    }
	return rf;
}




