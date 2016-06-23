#include "FFreePeak_OnePeak.h"

extern mutex g_mutexStream;
//NT=1 NR=1 ND=2 PN=FUN_FreePeak_OnePeak AN=ALG_PSO_G PS=10 ME=5000 NP=3 NB=3 NO=1 IT1=P1+1#P2+1#P4+1# IT2=H1+100#H2+50#H3+50# IT3=B1+0.1#B2+0.3#B3+0.6#

FFreePeak_OnePeak::FFreePeak_OnePeak(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), v[param_numObj]), FreePeak(v){

	if (v.find(param_numGOpt) != v.end()) m_numGOpt = v[param_numGOpt]<m_numBox ?(int) v[param_numGOpt] : m_numBox;
	else m_numGOpt = m_numBox;

	initialize();
	m_tag.insert(MMP);
	setOptType(MAX_OPT);
	setObjSet();

	m_proPar << "Number of Global Opt.: " << m_numGOpt << ";";

	//m_initialBox4Sol.push_back(1);
	if (m_simulateFun) 	simulateFun();
}

void FFreePeak_OnePeak::createPeaks_(){
	
	string name;
	if (IS_PROBLEM_NAME(m_id, "FUN_FreePeak_OnePeak")){
		name = "FUN_OnePeak";
	}
	else if (IS_PROBLEM_NAME(m_id, "FUN_FreePeak_C_OnePeak")){
		name = "FUN_C_OnePeak";
	}
	
	struct Peak{
		vector<double> loc;
		int shape;
		double hei, minH;
		vector<int> transf;
	};
	vector<Peak> peak(m_numBox);
	for (auto&i : peak) i.loc.resize(m_numDim);
	//Global::g_arg[param_dataFile3] = "location_2_dim2_peak20.loc";
	if (Global::g_arg.find(param_dataFile3) != Global::g_arg.end()){	
		stringstream ss;
		ss << Global::g_arg[param_workingDir] << Global::g_arg[param_dataDirectory1] << Global::g_arg[param_dataFile3];
		ifstream in(ss.str().c_str());
		if (!in) throw myException("file does not exist @ FFreePeak_OnePeak::createPeaks_()");
		int numPeak = 0, numDim = 0;
		in >> numPeak >> numDim;
		if (numPeak != m_numBox||numDim!=m_numDim) throw myException("data do not match @ FFreePeak_OnePeak::createPeaks_()");
		for (int i = 0; i < m_numBox; ++i){
			for (int j = 0; j < m_numDim; ++j){
				in >> peak[i].loc[j];
			}
		}
		in.close();
	}
	else throw myException("No data available @ createPeaks_()");

	//Global::g_arg[param_dataFile1] = "config_94_peak20.conf";
	if (Global::g_arg.find(param_dataFile1) != Global::g_arg.end()){
		set<int> remainBox;
		set<int> used;
		for (int i = 0; i < m_numBox; ++i) remainBox.insert(i);
		stringstream ss;
		ss << Global::g_arg[param_workingDir] << Global::g_arg[param_dataDirectory1] << Global::g_arg[param_dataFile1];

		ifstream in(ss.str().c_str());
		if (!in) throw myException("file does not exist @ FFreePeak_OnePeak::createPeaks_()");

		int numPeak = 0,idx=0,numObj;
		string s;
		getline(in, s);
		getline(in, s);
		istringstream iss(s);
		iss >> numPeak >> numObj;
		if (numPeak != m_numBox) throw myException("data do not match @ FFreePeak_OnePeak::createPeaks_()");
		
		while (!in.eof()){
			string s;
			getline(in, s);
			if (s == "#END") break;

			vector<string> col(5);
			stringstream ss(s);
			for (int i = 0; i < col.size(); ++i) ss >> col[i];
			vector<TypeVar> value;
			int bidx;
			if (col[0] == "smallest"){
				bidx = smallestBox(used);
			}
			else if (col[0] == "largest"){
				bidx = largestBox(used);
			}
			else if (col[0] == "random"){
				set<int>::iterator it = remainBox.begin();
				int i = Global::msp_global->getRandInt(0, remainBox.size(), Program_Problem);
				for (int j = 0; j < i; j++) ++it;
				bidx = *(it);
			}
			else{
				gGetStringValue(col[0], value);
				bidx = value[0];
			}
			remainBox.erase(bidx);
			used.insert(bidx);

			gGetStringValue(col[1], value);
			peak[bidx].shape = value[0];

			gGetStringValue(col[2], value);
			if (value[0].is_int())	peak[bidx].hei = (int)value[0];
			else peak[bidx].hei = value[0];

			gGetStringValue(col[3], value);
			if (value[0].is_int())	peak[bidx].minH = (int)value[0];
			else peak[bidx].minH = value[0];
			

			gGetStringValue(col[4], value);
			vector<int> transf(value.size());
			for (int i = 0; i < transf.size(); ++i) transf[i] = value[i];

			peak[bidx].transf.resize(4);
			for (auto&j : peak[bidx].transf) j = 0;
			peak[bidx].transf[0] = 1;
			for (auto &j : transf){
				if (j < 3) peak[bidx].transf[0] = j + 1;
				else peak[bidx].transf[j - 2] = 1;
			}
		}

		m_numPeak = m_numBox;
		int idxp = 0, idxb = 0;
		for (int t = 0; t < m_tree.size(); ++t){
			for (int i = 0; i < m_tree[t]->region.size(); ++i, ++idxb){
				(*m_box)[idxb].peak.resize(m_peaksPerBox);
				for (int j = 0; j < m_peaksPerBox; ++j){
					FOnePeak *p = nullptr;
					p = new FOnePeak(Global::msm_pro[name], m_numDim, name, peak[idxb].loc, peak[idxb].shape, peak[idxb].hei, &peak[idxb].transf, idxp);
					p->gidx() = i;
					p->idx() = j;
					p->tidx() = t;
					p->setBasinRatio((*m_box)[idxb].boxRatio);
					p->setMaxHeight(mc_maxHeight);
					p->setStdMinHeight(peak[idxb].minH);
					p->setSdandardize(true);
					(*m_box)[idxb].peak[j].reset(p);
					m_fp.push_back(p);
					m_objIdx[make_tuple(idxb, j, 0)] = 0;
					m_peakIdx[make_pair(idxb, j)] = idxp;
					++idxp;
				}
			}
		}
		
	}
	else throw myException("No data available @ createPeaks_()");
}

ReturnFlag FFreePeak_OnePeak::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode , bool flag ){
	CodeVReal &s = dynamic_cast< CodeVReal&>(ss);
	int tidx = getTreeRoot(s.m_x);
	int gidx = m_tree[tidx]->get_regionIdx(s.m_x);
	int bidx = boxIdx(tidx, gidx);
	for (size_t i = 0; i < (*m_box)[bidx].peak.size(); ++i){
		if (i>0 && ss.m_obj[0] >= dynamic_cast<FOnePeak*>((*m_box)[bidx].peak[i].get())->height()) continue;
		CodeVReal x(mapToPeak(s.m_x, gidx, i,tidx), s.m_obj);
		(*m_box)[bidx].peak[i]->evaluate_(x, rFlag, mode, false);

#ifdef OFEC_DEBUG_
		if (mode == Program_Algorithm&&Global::msp_global->mp_problem&&rFlag&&dynamic_cast<FOnePeak*>((*m_box)[bidx].peak[i].get())->height() - x.m_obj[0] <= m_accuracy) {
			(*m_box)[bidx].peak[i]->getGOpt().setFoundFlagTrue();
		}
#endif		
		if (i == 0) ss.m_obj[0] = x.m_obj[0];
		else if (ss.m_obj[0] < x.m_obj[0])		ss.m_obj[0] = x.m_obj[0];
	}
	if (rFlag)		m_evals++;
	if (mode == Program_Algorithm&&Global::msp_global->mp_problem&&!Global::msp_global->mp_problem->isProTag(MOP)) m_globalOpt.isFound(s, m_disAccuracy, m_accuracy);

	if (Global::msp_global->mp_algorithm != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()){ 
		//for debuging
		#ifdef OFEC_DEBUG_
		outPeakRate();
		#endif
		return Return_Terminate;
	}
	return Return_Normal;
}
void FFreePeak_OnePeak::setHeight(const vector<double>&h){
	double max = *max_element(h.begin(), h.end());
	for (int i = 0; i < m_numPeak; ++i){
		dynamic_cast<FOnePeak*>(m_fp[i])->setHeight(h[i]);
		//for test, warning, time consuming
#ifdef OFEC_DEBUG_
		dynamic_cast<FOnePeak*>(m_fp[i])->difficulty(h[i]/max);
#endif
	}
}
void FFreePeak_OnePeak::outPeakRate(){
	static vector<double> pr(m_numPeak,0);
	for (int i = 0; i < m_numPeak; ++i){
		if (dynamic_cast<FOnePeak*>(m_fp[i])->getGOpt().getNumGOptFound()==1) pr[i] += 1;
	}
	g_mutexStream.lock();
	int run = Global::g_arg[param_numRun];
	string ss = Global::g_arg[param_workingDir];
	ss += "Result/";
	ss += mSingleObj::getSingleObj()->m_fileName.str();
	ss += "PR.txt";
	ofstream out(ss.c_str());
	for (int i = 0; i < m_numPeak; ++i){
		out << pr[i] / run << " " << dynamic_cast<FOnePeak*>(m_fp[i])->mean() << " " << dynamic_cast<FOnePeak*>(m_fp[i])->variance()<<endl;
	}
	out.close();
	g_mutexStream.unlock();
	

}

void FFreePeak_OnePeak::simulateFun(){
	//simulate traditional functions using peak forest
	vector<double> h(m_numPeak);
	for (int i = 0; i < m_numPeak; ++i){
		const vector<double> & loc = m_peak[i].m_x;
		//set a peak height using a traditional function
		h[i] = 0;
		//the Sphere function
		/*for (int d = 0; d < loc.size(); ++d){
			h[i] += loc[d] * loc[d];
		}*/
		//the Schwefel function
		for (int d = 0; d < loc.size(); ++d)
		{
			double x = -500 + 1000 * (loc[d] + 100) / 200;
			h[i] += -x * sin(sqrt(fabs(x)));
		}
	}
	double maxh = *max_element(h.begin(), h.end());
	double minh = *min_element(h.begin(), h.end());
	for (int i = 0; i < m_numPeak; ++i){
		h[i] = 1 + 99 * (h[i] - minh) / (maxh - minh);
	}

	for (int i = 0; i < m_numPeak; ++i){
		dynamic_cast<FOnePeak*>(m_fp[i])->setHeight(h[i]);
	}

}
