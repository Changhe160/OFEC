#include "FFreePeak_M_OnePeak.h"
#include "../../FOnePeak.h"

mutex pf_mutex;
FFreePeak_M_OnePeak::FFreePeak_M_OnePeak(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), v[param_numObj]), FreePeak(v), \
m_numParetoRegion(v[param_numParetoRegion]), m_case(CASE((int)v[param_case])), m_pfr(m_numObj), m_isPR(m_numBox,false),m_attraction(true){
	m_peaksPerBox = m_numObj;
	//m_dm = FreePeak::DivisionMode::DM_Even;
	if (m_numParetoRegion > m_numBox) m_numParetoRegion = m_numBox;
	if (m_case == Jump || m_case == Web) {
		m_objShape.resize(m_numObj);
	}
	else {
		m_objShape.resize(m_numObj*m_numBox);
	}
	initialize();
	setOptType(MAX_OPT, -1);
	configAttraction();
	computePFR();
	computePS();	
	setObjSet();
}
bool FFreePeak_M_OnePeak::isValid(const VirtualEncoding  &s){
	if (m_name == "FUN_FreePeak_M_C_OnePeak2"){
		//a region within a sphere with radius r is valid 
		const CodeVReal  &x = dynamic_cast<const CodeVReal  &>(s);

		int tidx = getTreeRoot(x.m_x);
		int gidx = m_tree[tidx]->get_regionIdx(x.m_x);
		int bidx = boxIdx(tidx, gidx);

		bool flag = true;
		for (size_t i = 0; i < (*m_box)[bidx].peak.size(); ++i){
			flag = flag && (*m_box)[bidx].peak[i]->isValid(CodeVReal(mapToPeak(x.m_x, gidx, i,tidx),vector<double>(1)));
			if (!flag) return false;
		}
	}
	return true;
}

void FFreePeak_M_OnePeak::createPeaks_(){
	m_numPeak = 0;
	m_fp.clear();
	//m_objIdx.clear();
	m_peakIdx.clear();
	string name;
	if (m_name == "FUN_FreePeak_M_OnePeak"){
		name = "FUN_OnePeak";	
		//name = "FUN_R_OnePeak";	
	}
	else if (m_name == "FUN_FreePeak_M_C_OnePeak"){
		name = "FUN_C_OnePeak";
	}


	struct Peak{
		vector<double> loc;
		int shape;
		double hei, minH;
		vector<int> transf;
		double radius;
	};
	
	if (Global::g_arg.find(param_dataFile3) != Global::g_arg.end() && Global::g_arg.find(param_dataFile1) != Global::g_arg.end()){
		vector<Peak> peak;
		stringstream ss;
		ss << Global::g_arg[param_workingDir] << Global::g_arg[param_dataDirectory1] << Global::g_arg[param_dataFile3];
		ifstream in(ss.str().c_str());
		if (!in) throw myException("file does not exist @ FFreePeak_OnePeak::createPeaks_()");
		string str;
		getline(in, str);
		int numBox = 0, numDim = 0,numObj;
		getline(in, str);
		istringstream istr(str);		
		istr >> numBox >> numObj >> numDim;

		if (numBox != m_numBox || numDim != m_numDim||numObj!=m_numObj) throw myException("data do not match @ FFreePeak_OnePeak::createPeaks_()");
		peak.resize(m_numBox*m_numObj);
		for (auto&i : peak) i.loc.resize(m_numDim);
		getline(in, str);
		for (int i = 0; i < m_numBox*m_numObj; ++i){
			for (int j = 0; j < m_numDim; ++j){
				in >> peak[i].loc[j];
			}
			in >> peak[i].radius;
		}
		in.close();

		ss.str("");
		ss << Global::g_arg[param_workingDir] << Global::g_arg[param_dataDirectory1] << Global::g_arg[param_dataFile1];

		in.open(ss.str().c_str());
		if (!in) throw myException("file does not exist @ FFreePeak_OnePeak::createPeaks_()");

		getline(in, str);
		getline(in, str);
		istringstream iss(str);
		iss >> numBox >> numObj;
		vector<TypeVar> value;
		getline(in, str);
		getline(in, str);
		gGetStringValue(str, value);

		vector<int> paretobox;
		for (auto i : value) {
			paretobox.push_back(i);
			m_isPR[i] = true;
			m_ipsr1.push_back(i);
		}
		m_ipsr.push_back(paretobox);
		getline(in, str);
		if (numBox != m_numBox || numObj!=m_numObj) throw myException("data do not match @ FFreePeak_OnePeak::createPeaks_()");
		int pidx=0;
		size_t pos;
		for (int num = 0; num < m_numBox*m_numObj; ++num) {
			string s;
			getline(in, s);
			vector<string> col(5);
			stringstream ss(s);
			for (int i = 0; i < col.size(); ++i) ss >> col[i];
						
			gGetStringValue(col[1], value);
			
			peak[pidx].shape = value[0];
			if (m_case == Jump || m_case == Web) m_objShape[num%m_numObj] = value[0];
			else m_objShape[num] = value[0];

			gGetStringValue(col[2], value);
			if (value[0].is_int())	peak[pidx].hei = (int)value[0];
			else peak[pidx].hei = value[0];

			gGetStringValue(col[3], value);
			if (value[0].is_int())	peak[pidx].minH = (int)value[0];
			else peak[pidx].minH = value[0];

			gGetStringValue(col[4], value);
			vector<int> transf(value.size());
			for (int i = 0; i < transf.size(); ++i) transf[i] = value[i];

			peak[pidx].transf.resize(4);
			for (auto&j : peak[pidx].transf) j = 0;
			peak[pidx].transf[0] = 1;
			for (auto &j : transf){
				if (j < 3) peak[pidx].transf[0] = j + 1;
				else peak[pidx].transf[j - 2] = 1;
			}
			pidx++;
		}

		m_numPeak = 0;
		int bidx = 0;
		for (int t = 0; t < m_tree.size(); ++t){
			for (int i = 0; i < m_tree[t]->region.size(); ++i, ++bidx){
				(*m_box)[bidx].peak.resize(m_peaksPerBox);
				for (int j = 0; j < m_peaksPerBox; ++j){
					FOnePeak *p = 0;
					if (m_case == Jump || m_case == Web) {
						p = new FOnePeak(Global::msm_pro[name], m_numDim, name, peak[m_numPeak].loc, m_objShape[j], peak[m_numPeak].hei);
						if (j == 0 && (m_objShape[0] == FOnePeak::SH7 || m_objShape[0] == FOnePeak::SH8 || m_objShape[0] == FOnePeak::SH9 || m_objShape[0] == FOnePeak::SH10)) {
							p->setRadiusStepFun(peak[m_numPeak].radius);
						}
					}
					else p = new FOnePeak(Global::msm_pro[name], m_numDim, name, peak[m_numPeak].loc, m_objShape[m_numPeak], peak[m_numPeak].hei);
					
					p->gidx() = i;
					p->idx() = j;
					p->tidx() = t;
					p->setSdandardize(true);
					if (name == "FUN_C_OnePeak"){
						p->setValidRadius(Global::g_arg[param_validRadius]);
					}
					(*m_box)[bidx].peak[j].reset(p);
					m_fp.push_back(p);
					//m_objIdx[make_tuple(bidx, j, 0)] = m_numPeak;
					m_peakIdx[make_pair(bidx, j)] = m_numPeak;
					++m_numPeak;
				}
			}
		}

	}else{
		throw myException("No data available @ createPeaks_()");
	}
	
	
}

ReturnFlag FFreePeak_M_OnePeak::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag){
	CodeVReal &s = dynamic_cast< CodeVReal&>(ss);
	int tidx = getTreeRoot(s.m_x);
	int gidx = m_tree[tidx]->get_regionIdx(s.m_x);
	int bidx = boxIdx(tidx, gidx);

	for (size_t i = 0; i < (*m_box)[bidx].peak.size(); ++i){
		CodeVReal x(mapToPeak(s.m_x, gidx, i,tidx), 1);
		(*m_box)[bidx].peak[i]->evaluate_(x, rFlag, mode, false);
		ss.m_obj[i] = x.m_obj[0];
	}
	if (rFlag)		m_evals++;

	if (mode==Program_Algorithm&& Global::msp_global->mp_algorithm != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()){ return Return_Terminate; }
	return Return_Normal;
}

void FFreePeak_M_OnePeak::configAttraction(){
	
	switch (m_case)	{
	case Jump:
		if (!m_attraction){			
			for (int i = 0; i < m_numBox; ++i){
				if (m_isPR[i]) continue;
				for (int j = 1; j < m_peaksPerBox; ++j){
					FOnePeak *p = dynamic_cast<FOnePeak*>((*m_box)[i].peak[j].get());
					double h = p->height();
					p->setHeight(h*0.99);
				}
			}
		}
		break;
	case Web:
		if (m_attraction){
			for (int i = 0; i < m_numBox; ++i){
				int so = Global::msp_global->getRandInt(1, m_peaksPerBox);
				for (int j = 1; j < m_peaksPerBox; ++j){
					if (j != so){
						vector<double> &loc = dynamic_cast<FOnePeak*>((*m_box)[m_ipsr[0][0]].peak[j].get())->location();
						FOnePeak *p = dynamic_cast<FOnePeak*>((*m_box)[i].peak[j].get());
						p->setLocation(loc);
					}
				}
			}
		}		
		break;
	default:	break;
	}
}
void FFreePeak_M_OnePeak::computePFR(){

	switch (m_case)
	{
	case Jump:	case Countable:	case Web:
		for (int i = 0; i < m_ipsr.size(); ++i){
			vector<vector<POS>> ps;
			for (int k = 0; k < m_ipsr[i].size(); ++k){
				vector<POS> psr;
				for (int j = 0; j < m_numObj; ++j){
					int idx = m_peakIdx[make_pair(m_ipsr[i][k], j)];
					psr.push_back(m_peak[idx]);
					psr.back().isboundary = true;
					psr.back().ispeak = true;
				}
				ps.push_back(move(psr));
			}
			m_psr.push_back(move((ps)));
		}
		break;
	default:
		break;
	}
	

}
void FFreePeak_M_OnePeak::getPSBetween(const POS &p1, POS &p2, list<POS>& ps, FFreePeak_M_OnePeak&p, int maxNum, int oidx, int order, int region, bool ispeak, bool included){
	
	double precesion = 1.E-5;
	
	double dpf = 0.5;  // density
	list<POS> nps;
	bool good = false;
	int increase = -1;
	
	do{
		nps.clear();
		int num = order;
		nps.push_back(p1);
		nps.back().order = num;
		nps.back().isboundary = ispeak;
		nps.back().pfr = region;
		bool finish = false, nopoint = true;
		CodeVReal x1 = p1.sol, x2 = p2.sol, x3, x4(p.m_numDim, p.m_numObj);
		do{			
			nopoint = true;
			while (x2.getObjDistance(x1.m_obj, DIS_EUCLIDEAN) >= dpf){
				x3 = x2;
				for (int i = 0; i < p.m_numDim; ++i){
					x2.m_x[i] = (x1.m_x[i] + x2.m_x[i]) / 2;
				}
				p.evaluate_(x2, false, Program_Problem, false);
				nopoint = false;
			}
			if (!nopoint){
				double dis;
				do{
					for (int i = 0; i < p.m_numDim; ++i){
						x4.m_x[i] = (x3.m_x[i] + x2.m_x[i]) / 2;
					}
					p.evaluate_(x4, false, Program_Problem, false);
					dis = x4.getObjDistance(x1.m_obj, DIS_EUCLIDEAN);
					if (dis > dpf){
						x3 = x4;
					}
					else{
						x2 = x4;
					}
				} while (fabs(dis - dpf) >= precesion);

				++num;
				nps.push_back(x4);
				nps.back().isboundary = ispeak;
				nps.back().order = num;
				nps.back().pfr = region;
				x1 = x4;
				x2 = p2.sol;
			}
			else{
				/*for (int i = 0; i < p.m_numDim; ++i){
					x4.m_x[i] = (x1.m_x[i] + x2.m_x[i]) / 2;
				}
				p.evaluate_(x4, false, Program_Problem, false);*/
				finish = true;
			}

		} while (!finish);
		++num;
		if (!included){
			nps.push_back(p2);
			nps.back().isboundary = ispeak;
			nps.back().order = num;
			nps.back().pfr = region;
		}
		else{
			p2.order = num;
			p2.pfr = region;
		}
		if (!good&&nps.size() > maxNum){
			dpf *= 2;
			increase = 0;
		}
		else if (increase == 0 || nps.size() == maxNum){ good = true; }
		
		if (!good&&nps.size() < maxNum){
			dpf /= 2;
			increase = 1;
		}
		else if (increase == 1 || nps.size() == maxNum){ good = true; }

	} while (!good);
	
	
	ps.splice(ps.end(), nps);
}


void FFreePeak_M_OnePeak::computePS(){

	m_ps.clear();
	if (m_case != Countable) {
		for (int i = 0; i < m_psr.size(); ++i) {
			list<POS> psa;
			for (int j = 0; j < m_psr[i].size(); ++j) {
				list<POS> ps, hdps;
				int objs = m_psr[i][j].size();

				POS x(m_psr[i][j][objs - 1]);
				x.ispeak = true;
				x.isboundary = true;
				ps.push_back(move(x));

				for (int k = objs - 2; k >= 0; --k) {
					list<POS>::iterator first = --ps.end();
					int num = ps.size(), n = 0;

					int numTask = std::thread::hardware_concurrency();
					if (numTask>num) numTask = num;

					vector<thread> atrd;
					vector<future<int>> af(numTask);

					size_t unit = (num / numTask);
					vector<vector<list<POS>::iterator>> tsk(numTask);
					for (auto &t : tsk)	t.reserve(unit + 1);

					size_t ii = 0;
					list<POS>::iterator it = ps.begin();
					for (; ii < num; ++ii, ++it) {
						int k = ii / unit;
						if (k == numTask) break;
						tsk[k].push_back(it);
					}
					for (size_t jj = 0; ii < num; ++ii, ++jj, ++it) {
						tsk[jj].push_back(it);
					}
					for (ii = 0; ii<numTask; ii++) {
						FOnePeak*peak = dynamic_cast<FOnePeak*>((*m_box)[m_ipsr[i][j]].peak[k].get());
						packaged_task<int(int oidx, int shape, POS &o, list<POS>& ps, vector<list<POS>::iterator> &tsk, FFreePeak_M_OnePeak&p)> at(PSThread);
						af[ii] = at.get_future();
						atrd.push_back(thread(move(at), k, peak->shape(), std::ref(m_psr[i][j][k]), std::ref(ps), std::ref(tsk[ii]), std::ref(*this)));
					}
					for (auto&t : atrd) t.join();
					for (auto&f : af) f.wait();

					first++;
					first->ispeak = true;
					first->isboundary = true;


					//make the POF uniformally distributed
					vector<vector<list<POS>::iterator>> sameOrder;
					for (auto it = ps.begin(); it != ps.end(); ++it) {
						if (it->order >= sameOrder.size()) {
							sameOrder.resize(it->order + 1);
						}
						sameOrder[it->order].push_back(it);
					}
					int maxItem = 0;
					int maxIdx = 0;
					for (int it = 0; it< sameOrder.size(); ++it) {
						if (maxItem <= sameOrder[it].size()) {
							maxItem = sameOrder[it].size();
							maxIdx = it;
						}
					}

					vector<list<POS>::iterator> stay;
					double step = maxItem*1.0 / (maxIdx + 1);
					for (auto it = 0; it< maxIdx; ++it) {
						int num = ceil((it)*step);
						vector<list<POS>::iterator> remain, keep;

						if (it == 0) {
							for (auto init = sameOrder[it].begin(); init != sameOrder[it].end(); ++init) {
								if ((*init)->ispeak) {
									keep.push_back(*init);
									break;
								}
							}

						}
						else {
							for (auto init = sameOrder[it].begin(); init != sameOrder[it].end(); ++init) {
								if ((*init)->isboundary&&!(*init)->ispeak) {
									keep.push_back(*init);
								}
								else  remain.push_back(*init);
							}


							if (keep.size() < num) {
								if (remain.size() > num) {
									vector<int> ran;
									vector<int> total(remain.size());
									for (int ii = 0; ii < total.size(); ++ii) total[ii] = ii;
									for (int ii = 0; ii < num - keep.size(); ++ii) {
										int r = Global::msp_global->mp_uniformPro->Next()*total.size();
										ran.push_back(total[r]);
										total.erase(total.begin() + r);
									}
									for (auto ii : ran) {
										keep.push_back(remain[ii]);
									}
								}
								else {
									for (auto ii : remain) {
										keep.push_back(ii);
									}
								}
							}
						}

						for (auto &ii : keep) stay.push_back(ii);
					}

					// put remaining points in ps
					for (auto it = maxIdx; it < sameOrder.size(); ++it) {
						for (auto &jj : sameOrder[it]) {
							stay.push_back(jj);
						}
					}
					for (auto &ii : stay) {
						psa.push_back(*ii);
					}
				}
				//psa.splice(psa.end(), ps);
			}
			m_ps.splice(m_ps.end(), psa);
		}
	}
	else {
		for (int b = 0; b < m_numBox; ++b) {
			if (m_isPR[b]) {
				m_ps.push_back(POS(m_peak[m_peakIdx[make_pair(b,0)]]));
			}
		}
	}	
	m_globalOpt.setFlagLocTrue();
	m_globalOpt.setNumOpts(m_ps.size());
	int j = 0;
	for (auto i = m_ps.begin(); i != m_ps.end();++i,++j){
		m_globalOpt[j] = i->sol;
	}
}
void FFreePeak_M_OnePeak::getSegmentPS(int oidx, int shape, const POS &o1, const POS &o2, vector<pair<POS, POS>> &mseg){
	mseg.clear();
	vector<pair<vector<double>, vector<double>>> seg;
	if (shape == FOnePeak::SH8||shape==FOnePeak::SH9){	
		int gidx = m_tree[0]->get_regionIdx(o1.sol.m_x);		
		FOnePeak*peak = dynamic_cast<FOnePeak*>((*m_box)[gidx].peak[oidx].get());
		peak->getContiSegment(mapToPeak(o2.sol.m_x, gidx, oidx), seg);
		for (auto &s : seg){
			CodeVReal p1(m_numDim, m_numObj), p2(m_numDim, m_numObj);
			p1.m_x = mapFromPeak(s.first, gidx, oidx);
			p2.m_x = mapFromPeak(s.second, gidx, oidx);
			evaluate_(p1, false, Program_Problem, false);
			evaluate_(p2, false, Program_Problem, false);
			mseg.push_back(make_pair(POS(p1), POS(p2)));
		}
	
	}
	else if (m_objShape[0] == FOnePeak::SH8){
		int gidx = m_tree[0]->get_regionIdx(o1.sol.m_x);
		FOnePeak*peak = dynamic_cast<FOnePeak*>((*m_box)[gidx].peak[0].get());
		peak->getContiSegment(mapToPeak(o1.sol.m_x, gidx, 0),mapToPeak(o2.sol.m_x, gidx, 0), seg);
		for (auto &s : seg){
			CodeVReal p1(m_numDim, m_numObj), p2(m_numDim, m_numObj);
			p1.m_x = mapFromPeak(s.first, gidx, 0);
			p2.m_x = mapFromPeak(s.second, gidx, 0);
			evaluate_(p1, false, Program_Problem, false);
			evaluate_(p2, false, Program_Problem, false);
			mseg.push_back(make_pair(POS(p1), POS(p2)));
		}
	}
	else{
		mseg.push_back(make_pair(o1, o2));
	}


}
int FFreePeak_M_OnePeak::PSThread(int oidx, int shape, POS &o, list<POS>& ps, vector<list<POS>::iterator> &tsk, FFreePeak_M_OnePeak&p){

	for (auto t : tsk){
		list<POS> hdps;
	
		vector<pair<POS, POS>> mseg;
		p.getSegmentPS(oidx, shape, o, *t, mseg);
		

		vector<int> maxNum(p.m_numObj), initNum = { 50, 25, 10, 3};
		if (p.m_numObj <= 3) {
			initNum[0]=initNum[1] = 200;
		}
		for (int i = 0; i < p.m_numObj; ++i){
			if (i < initNum.size()) maxNum[i] = initNum[i];
			else maxNum[i] = 3;
		}
		int curMaxNum = maxNum[p.m_numObj - oidx - 2];
		vector<int> segMaxNum(mseg.size());
		
		if (p.m_objShape[0] == 8){
			int remain = curMaxNum;
			for (auto i = 0; i < segMaxNum.size(); ++i){
				if (i == segMaxNum.size() - 1) segMaxNum[i] = remain;
				else	segMaxNum[i] = remain/2;
				remain -= segMaxNum[i];
				if (segMaxNum[i]  < 3) segMaxNum[i] = 3;
			}
		}
		else{
			for (auto i = 0; i < segMaxNum.size(); ++i){
				segMaxNum [i]= ceil(curMaxNum*1.0 / mseg.size());
				if (segMaxNum[i]  < 3) segMaxNum[i] = 3;
			}
			
		}
		

		int i = 0,region=0;
		for (auto s = mseg.begin(); s != mseg.end(); ++s, ++region){
			if (s == mseg.end() - 1)	getPSBetween(s->first, s->second, hdps, p, segMaxNum[region], oidx, i, region, t->ispeak, true);
			else getPSBetween(s->first, s->second, hdps, p, segMaxNum[region], oidx, i, region, t->ispeak, false);
			i += hdps.size();
			{
				unique_lock<mutex> lock(pf_mutex);
				ps.splice(ps.end(), hdps);
			}
		}
		t->order = i;
		t->pfr = region - 1;
	}
	return 0;
}

const vector<int>& FFreePeak_M_OnePeak::getParetoRegion(){
	return m_ipsr1;
}