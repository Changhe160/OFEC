#include "ZDT.h"

ZDT::ZDT(ParamMap &v):BenchmarkFunction((v[param_proId]),(v[param_numDim]),(v[param_proName]),2)
{
    setSearchRange(0.,1.);
	vector<ProTag> p_tag(1,MOP);
	p_tag.push_back(CONT);
	setProTag(p_tag);
	setOptType(MIN_OPT,-1);
	m_popInitialMode=POP_INIT_UNIFORM;
}

void ZDT::generateAdLoadPF()
{
	int num=1000;
	m_globalOpt.setNumOpts(num);
	m_originalGlobalOpt.setNumOpts(num);
	m_originalGlobalOpt.setFlagLocTrue();

	ofstream out;
	ifstream infile;
	stringstream os;
	os<<Global::msp_global->g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/PF_"<<Global::g_arg[param_proName]<<"_Opt_"<<num<<".txt";
	infile.open(os.str());
	if(!infile)
	{
		out.open(os.str());
		double temp;
		for(int i=0;i<num;i++)
		{
			temp=static_cast<double>(i)/num;
			m_originalGlobalOpt[i].data().m_x[0]=temp;
			for(int j=1;j<m_numDim;j++)
				m_originalGlobalOpt[i].data().m_x[j]=0.;
			BenchmarkFunction::evaluate_(m_originalGlobalOpt[i].data(),false);
			out<<m_originalGlobalOpt[i].obj(0)<<" "<<m_originalGlobalOpt[i].obj(1)<<endl;
		}
		out.close();
		m_globalOpt=m_originalGlobalOpt;
	}
	else 
	{
		double temp;
		for(int i=0;i<num;i++)
		{
			temp=static_cast<double>(i)/num;
			m_originalGlobalOpt[i].data().m_x[0]=temp;
			for(int j=1;j<m_numDim;j++)
				m_originalGlobalOpt[i].data().m_x[j]=0.;
			infile>>m_originalGlobalOpt[i].obj()[0];
			infile>>m_originalGlobalOpt[i].obj()[1];
		}
		m_globalOpt=m_originalGlobalOpt;
		infile.close();
	}
}
