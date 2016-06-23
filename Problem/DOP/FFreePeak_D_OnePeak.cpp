#include "FFreePeak_D_OnePeak.h"
#include "../FunctionOpt/FOnePeak.h"
#ifdef OFEC_DEMON
#include "../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
#endif
FFreePeak_D_OnePeak::FFreePeak_D_OnePeak(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1), \
FreePeak(v), m_feature1(static_cast<Feature1>((int)v[param_changeType])), m_noise(v[param_noiseFlag]), m_numPeakChange(v[param_flagNumPeakChange]),\
m_feature2(static_cast<Feature2>((int)v[param_changeType2])), m_shiftSeverity(v[param_shiftLength]), m_changeInterval(v[param_changeFre]),\
m_detect(v[param_changeRatio]){
	
	if (m_numPeakChange){
		m_peaksPerBox = v[param_peaksPerBox] = 1;
		m_timewindow = v[param_timeWindow] = 0;
	}
	else {
		if(v.find(param_timeWindow)!=v.end())		m_timewindow = v[param_timeWindow];
		else v[param_timeWindow]=m_timewindow;
	}
		
	setShapeSet();
	initialize();

	setAccuracy(1.e-2);
	m_tag.insert(DOP);
	m_tag.insert(MMP);
	m_initNumPeak = m_numPeak;

	if (m_feature1 == FT_default){
		m_changingPeak.resize(m_numPeak);
		for (int i = 0; i < m_numPeak; ++i) m_changingPeak[i] = i;
	}
	if (m_feature2 == FT_basin) uf1.m_basinChange = true;
	else{ 
		uf1.m_predic = v[param_predicFlag]; 
		if (uf1.m_predic) m_lambada = 1;
	}

	vector<string> f1 = { "","timelink","detect","default" }, f2 = { "basin", "predic", "default" };
	m_proPar << "Noise: " << m_noise << "; Number of Peak change: " << m_numPeakChange << "; feature1: " << f1[m_feature1]\
		<< "; feature2: " << f1[m_feature2] << "; Shift Severity: " << m_shiftSeverity << "; Change Interval: " << m_changeInterval
		<< "; Change ratio: " << m_detect << ";";
}

void FFreePeak_D_OnePeak::createPeaks_(){
	m_numPeak = 0;
	m_fp.clear();
	//m_objIdx.clear();
	m_peakIdx.clear();
	string name;
	if (IS_PROBLEM_NAME(m_id, "FUN_FreePeak_D_OnePeak")){
		name = "FUN_OnePeak";
	}
	else if (IS_PROBLEM_NAME(m_id, "FUN_FreePeak_D_R_OnePeak")){
		name = "FUN_R_OnePeak";
	}
	else if (IS_PROBLEM_NAME(m_id, "FUN_FreePeak_D_C_OnePeak")){
		name = "FUN_C_OnePeak";
	}
	for (int t = 0, bidx = 0; t < m_tree.size(); ++t){
		for (int i = 0; i < m_tree[t]->region.size(); ++i, ++bidx){
			int num = Global::msp_global->getRandInt(1, m_peaksPerBox + 1, Program_Problem);
			if (num>m_peaksPerBox) num = m_peaksPerBox;
			(*m_box)[bidx].peak.resize(num);
			for (int j = 0; j < num; ++j){
				FOnePeak *p = nullptr;
				int shape = m_shapeSet[Global::msp_global->getRandInt(0, m_shapeSet.size(), Program_Problem) % m_shapeSet.size()];
				if (num == 1)			p = new FOnePeak(Global::msm_pro[name], m_numDim, name,false,shape);
				else p = new FOnePeak(Global::msm_pro[name], m_numDim, name, true,shape);
				p->setSdandardize(true);
				p->setShiftSeverity(m_shiftSeverity);
				p->setHeightSeverity(mc_maxHeightSeverity*Global::msp_global->mp_uniformPro->Next());
				p->gidx() = i;
				p->idx() = j;
				p->tidx() = t;
				p->setMaxHeight(mc_maxHeight);
				p->setMemorySize(m_timewindow);
				(*m_box)[bidx].peak[j].reset(p);
				m_fp.push_back(p);
				//m_objIdx[make_tuple(bidx, j, 0)] = 0;
				m_peakIdx[make_pair(bidx, j)] = m_numPeak;
				++m_numPeak;
			}
		}
	}

	m_isFound.resize(m_numPeak, false);
	m_tl.resize(m_numPeak, false);

	setupGOpt();
}
void FFreePeak_D_OnePeak::updateChangingPeaksIdx(int &numChangeGOpt){
	numChangeGOpt = 0;
	if (m_feature1 == FT_timelink){
		for (int i = 0; i < m_numPeak; ++i)	{
			if (m_tl[i])	
				m_changingPeak.push_back(i);
		}
	}
	else if (m_feature1 == FT_detect){
		vector<int> ridx(m_numPeak);
		for (int i = 0; i < m_numPeak; ++i) ridx[i] = i;
		int num = static_cast<int>(m_detect*m_numPeak) > 0 ? static_cast<int>(m_detect*m_numPeak) : 1;
		for (int i = 0; i < num; ++i){
			int idx = static_cast<int>(Global::msp_global->mp_uniformPro->Next()*ridx.size());
			m_changingPeak.push_back(ridx[idx]);
			ridx.erase(ridx.begin() + idx);
		}
	}

	if (m_feature1 == FT_timelink || m_feature1 == FT_detect){
		for (size_t i = 0; i < m_changingPeak.size(); ++i){
			if (m_peak[m_changingPeak[i]].m_obj[0] == m_maxHeight) ++numChangeGOpt;
		}
	}
	else{
		numChangeGOpt = m_numGOpt;
	}
}

void FFreePeak_D_OnePeak::updateChangingBasinIdx(int &numChangeGOpt){
	numChangeGOpt = 0;
	if (m_feature1 == FT_timelink){
		for (int i = 0; i < m_numPeak; ++i)	{
			if (m_tl[i]){
				m_changingPeak.push_back(i);
				changeBasin(i);
			}
		}
	}
	else if (m_feature1 == FT_detect){
		vector<int> ridx(m_numPeak);
		for (int i = 0; i < m_numPeak; ++i) ridx[i] = i;
		int num = static_cast<int>(m_detect*m_numPeak) > 0 ? static_cast<int>(m_detect*m_numPeak) : 1;
		for (int i = 0; i < num; ++i){
			int idx = static_cast<int>(Global::msp_global->mp_uniformPro->Next()*ridx.size());
			m_changingPeak.push_back(ridx[idx]);
			changeBasin(ridx[idx]);
			ridx.erase(ridx.begin() + idx);
		}
	}

	if (m_feature1 == FT_timelink || m_feature1 == FT_detect){
		for (size_t i = 0; i < m_changingPeak.size(); ++i){
			if (m_peak[m_changingPeak[i]].m_obj[0] == m_maxHeight) ++numChangeGOpt;
		}
	}
	else{
		numChangeGOpt = m_numGOpt;
		for (int i = 0; i < m_numBox; ++i)	{
			changeBasin(i);			
		}
	}
}
void FFreePeak_D_OnePeak::changeHeight(int &numChangeGOpt){

	for (size_t i = 0; i < m_changingPeak.size(); ){
		FOnePeak*p = dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[i]]);
		double h = p->height();
		double delta = p->getHeightSeverity()*Global::msp_global->mp_normalPro->Next();
		if (h + delta < 0 || h + delta >= mc_maxHeight) h -= delta;
		else h += delta;

		p->changeHeight(h);
		++i;
	}
	updateMaxHeight();

}

void FFreePeak_D_OnePeak::changeLocation(){
	if (m_feature2==FT_basin){
		set<int> ts;
		for (auto it=m_basin.begin(); it != m_basin.end(); ++it){
			double mag = (it->high - it->low)*m_shiftSeverity / 2;
			double delta = Global::msp_global->getRandFloat(-mag, mag, Program_Problem);
			int tidx,gidx;
			treeIdx(it->idx, tidx, gidx);
			ts.insert(tidx);
			double val = m_divisionPoint[tidx][gidx][it->cutdim];
			if (val + delta>it->high || val + delta < it->low){
				val -= delta;
			}
			else{
				val += delta;
			}
			m_divisionPoint[tidx][gidx][it->cutdim] = val;
		}
		for (auto i : ts){
			m_tree[i]->buildIndex();
		}
		computeBoxSize();
		setDisAccuracy(disNearestPeak() / 2);
	}
	else{
		MyVector v(m_numDim);
		for (size_t p = 0; p < m_changingPeak.size(); ++p){
			v.zero();
			if (dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[p]])->preLocSize() > 0){
				for (int i = 0; i < m_numDim; ++i){
					v[i] = dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[p]])->location(i) - dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[p]])->preLoc()[i];
				}
			}
			if (v.length() == 0) v.norRandomize(Program_Problem);
			v.normalize();
			MyVector r(m_numDim);
			r.norRandomize(Program_Problem);

			for (int i = 0; i < m_numDim; ++i){
				v[i] = m_lambada*v[i] + (1 - m_lambada)*r[i];
			}
			v.normalize();
			v *= m_fp[m_changingPeak[p]]->getSearchRange().getDomainSize()*m_shiftSeverity;
			v += dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[p]])->location();
			dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[p]])->shiftLocation(v.data());
		}
	}
	setDisAccuracy(disNearestPeak() / 2);
}

void FFreePeak_D_OnePeak::addNoise(vector<double>& x, int boxIdx,int tidx,int gidx){
	double x_,r;
	r = (*m_box)[boxIdx].boxRatio;
	for (int i = 0; i < m_numDim; ++i){
		x_ = m_noiseSeverity*r*Global::msp_global->mp_uniformPro->Next();
		if (x_ < m_tree[tidx]->region[gidx].box[i].first) x_ = m_tree[tidx]->region[gidx].box[i].first;
		if (x_ > m_tree[tidx]->region[gidx].box[i].second) x_ = m_tree[tidx]->region[gidx].box[i].second;
		x[i] = x_;
	}
}

ReturnFlag FFreePeak_D_OnePeak::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag_){
	
	CodeVReal &s = dynamic_cast<CodeVReal&>(ss);
	
	int tidx = getTreeRoot(s.m_x);
	int gidx = m_tree[tidx]->get_regionIdx(s.m_x);
	int bidx = boxIdx(tidx, gidx);

	if (m_noise){
		addNoise(s.m_x, bidx,tidx,gidx);
	}

	for (size_t i = 0; i < (*m_box)[bidx].peak.size(); ++i){
		if (i>0 && ss.m_obj[0] >= dynamic_cast<FOnePeak*>((*m_box)[bidx].peak[i].get())->height()) continue;
		CodeVReal x (mapToPeak(s.m_x, gidx, i,tidx),s.m_obj);
		(*m_box)[bidx].peak[i]->evaluate_(x, rFlag, mode, false);
		if (i == 0) ss.m_obj[0] = x.m_obj[0];
		else if (ss.m_obj[0] < x.m_obj[0])		ss.m_obj[0] = x.m_obj[0];
		if (mode == Program_Algorithm&&rFlag){
			pair<int, int> pidx = make_pair(bidx, i);
			double dis = getDistance(s, m_peak[m_peakIdx[pidx]], DIS_EUCLIDEAN);
			if (s.getObjDistance(m_peak[m_peakIdx[pidx]].m_obj) <= m_accuracy&&dis<m_disAccuracy) m_isFound[m_peakIdx[pidx]] = true;
			if (dis < 0.01*(*m_box)[bidx].boxSize) m_tl[m_peakIdx[pidx]] = true;

			int nearest = nearestPeak(s);
			if (m_nearSolut.find(nearest) != m_nearSolut.end()){
				if (ss.m_obj[0] > m_nearSolut[nearest])	m_nearSolut[nearest] = ss.m_obj[0];
			}
			else{
				m_nearSolut[nearest] = ss.m_obj[0];
			}
		}		
	}

	if (mode == Program_Algorithm&&Global::msp_global->mp_problem&&!Global::msp_global->mp_problem->isProTag(MOP)) m_globalOpt.isFound(s, m_disAccuracy, m_accuracy);

	if (rFlag&&m_evals%m_changeInterval == 0){
		Solution<CodeVReal>::initilizeWB(s);
	}

	if (rFlag)    m_evals++;
	bool flag=true;
	#ifdef OFEC_CONSOLE
		if (Global::msp_global->mp_algorithm != nullptr)	flag = !Global::msp_global->mp_algorithm->ifTerminating();
		else flag = true;
	#endif
	
	if (rFlag&&m_evals%m_changeInterval == 0 && flag){	
		change();
		for (int p = 0; p < m_numPeak; ++p) {
			m_isFound[p] = false;
			m_tl[p] = false;
		}
		m_globalOpt.setFoundFlagFalse();
		if (m_feature2==FT_basin)		m_basin.clear();
		if (m_feature1 == FT_timelink || m_feature1 == FT_detect) m_changingPeak.clear();
		else{
			m_changingPeak.resize(m_numPeak);
			for (int i = 0; i < m_numPeak; ++i) m_changingPeak[i] = i;
		}
		if (mSingleObj::getSingleObj() != nullptr){			
			mSingleObj::getSingleObj()->addGOpt(Global::msp_global->m_runId, m_maxHeight);		
		}

	#ifdef OFEC_DEMON
			msp_buffer->updateFitnessLandsacpe_();
	#endif
	}
	ReturnFlag rf = Return_Normal;
	if (rFlag){
		if (Global::msp_global->mp_algorithm != nullptr){
			if (Global::msp_global->mp_algorithm->ifTerminating()){ rf = Return_Terminate; }
			else if (Global::msp_global->mp_problem->isProTag(DOP)){
				if ((m_evals + 1) % (m_changeInterval) == 0){
					rf = Return_ChangeNextEval;
				}
				if (m_evals % m_changeInterval == 0){
					rf = Return_Change;
				}
			}
		}
	}
	return rf;
}

void FFreePeak_D_OnePeak::change(){
	if (m_numPeakChange){
		changeNumBoxes();
	}
	else{
		int numChangeGOpt;
		if (m_feature2 == FT_basin) updateChangingBasinIdx(numChangeGOpt);
		else updateChangingPeaksIdx(numChangeGOpt);
		changeShape();
		changeLocation();
		changeHeight(numChangeGOpt);		
		updatePeak();
	}	
	++m_counter;
}

void FFreePeak_D_OnePeak::setupGOpt(){
	
	for (int i = 0; i < m_numPeak; ++i){
		dynamic_cast<FOnePeak*>(m_fp[i])->setHeight(mc_maxHeight*Global::msp_global->mp_uniformPro->Next());
	}
	updateMaxHeight();
}

void FFreePeak_D_OnePeak::changeNumBoxes(){
	int newboxes;
	int period = 25, offset = m_initNumPeak / m_step;
	if ((m_counter + offset - m_numGOpt / m_step) / period % 2 == 0) newboxes = ((m_counter + offset - m_numGOpt / m_step) % period)*m_step;
	else newboxes = period*m_step - ((m_counter + offset - m_numGOpt / m_step) % 25)*m_step;

	newboxes += m_numGOpt; //minimum number of peaks is the numGOpt
	set<int> ts;
	//cout << "newbox: " << newboxes << "box: " << m_numBox<<endl;
	if (m_numBox > newboxes){ // randomly remove m_numPeak-newboxes peaks
		for (int i = 0; i < m_numBox - newboxes; ++i){
			int bidx = static_cast<int>(Global::msp_global->mp_uniformPro->Next()*(m_numBox-i));
			int tidx, gidx;
			treeIdx(bidx, tidx, gidx);
			ts.insert(tidx);
			m_divisionPoint[tidx].erase(m_divisionPoint[tidx].begin() + gidx);
		}
	}
	else if (m_numBox < newboxes){
		for (int i = 0; i <newboxes - m_numBox; ++i){
			int tidx;
			if (m_tree.size() == 1) tidx = 0;
			else tidx = Global::msp_global->getRandInt(0,m_tree.size(),Program_Problem);
			ts.insert(tidx);
			vector<double> p(m_numDim);
			for (auto i = 0; i < p.size();++i){
				p[i] = Global::msp_global->getRandFloat(m_tree[tidx]->get_rootBox().box[i].first, m_tree[tidx]->get_rootBox().box[i].second, Program_Problem);
			}
			m_divisionPoint[tidx].push_back(move(p));
		}
	}
	else return;
	
	m_numBox = newboxes;
	m_box.reset(new vector<Box>(m_numBox));
	for (auto i : ts){
		KDTreeSpace::PartitioningKDTree<double> *newtree = new KDTreeSpace::PartitioningKDTree<double>(m_numDim, m_divisionPoint[i]);
		newtree->setInitBox(m_tree[i]->get_rootBox().box);
		newtree->buildIndex();
		m_tree[i].reset(newtree);
	}
	computeBoxSize();	
	createPeaks();
}
void FFreePeak_D_OnePeak::changeBasin(const int bidx){
		
	BI b;
	int tidx, gidx;
	treeIdx(bidx,tidx,gidx);
	m_tree[tidx]->get_leafParent(gidx, b.idx, b.cutdim, b.low, b.high);
	vector<BI>::iterator it = m_basin.begin();
	for (; it != m_basin.end(); ++it){
		if (it->idx == b.idx) break; // already in m_basin
	}
	if (it == m_basin.end()) m_basin.push_back(move(b));
		
}

bool FFreePeak_D_OnePeak::predictChange(const int evalsMore){
	int evals = m_evals % m_changeInterval;
	if (evals + evalsMore >= m_changeInterval) return true;
	else return false;
}
double FFreePeak_D_OnePeak::getErr(){
	return m_maxHeight- Solution<CodeVReal>::getBestSolutionSoFar().data().m_obj[0];
}
double FFreePeak_D_OnePeak::getRobustErr(){
	double min=-1;
	for (int i = 0; i < m_numPeak; ++i){
		if (m_isFound[i]){
			double robust = dynamic_cast<FOnePeak*>(m_fp[i])->robustness();
			if (min == -1) min = robust;
			else{
				if (min>robust)min = robust;
			}
		}
	}
	return min;
}
int FFreePeak_D_OnePeak::nearestPeak(const CodeVReal &s){
	double dis = getDistance(m_peak[0], s,DIS_EUCLIDEAN);
	int idx = 0;
	for (int i = 1; i < m_numPeak; ++i){
		double d = getDistance(m_peak[i], s, DIS_EUCLIDEAN);
		if (dis > d){
			dis = d;
			idx = i;
		}
	}
	return idx;
}
double FFreePeak_D_OnePeak::getGPR(){
	int num = 0;
	for (int i = 0; i < m_numPeak; ++i){
		if (m_isFound[i]){
			if (dynamic_cast<FOnePeak*>(m_fp[i])->height() == m_maxHeight) num++;
		}
	}
	return num *1. / m_numGOpt;
}
double FFreePeak_D_OnePeak::getPR(){
	
	int num = 0;
	for (int i = 0; i < m_numPeak; ++i){
		if (m_isFound[i]){
			num++;
		}
	}
	return num*1. / m_numPeak;
}

int FFreePeak_D_OnePeak::getPT(){

	int num = 0;
	for (int i = 0; i < m_numPeak; ++i){
		if (m_isFound[i]){
			num++;
		}
	}
	return num;
}

bool FFreePeak_D_OnePeak::isAllGOptTraced(){
	int num = 0;
	for (int i = 0; i < m_numPeak; ++i){
		if (m_isFound[i]){
			if (dynamic_cast<FOnePeak*>(m_fp[i])->height() == m_maxHeight) num++;
		}
	}
	return (num==m_numGOpt)?true:false;
}

void FFreePeak_D_OnePeak::updateMaxHeight(){
	m_maxHeight = numeric_limits<double>::min();
	for (int i = 0; i < m_numPeak; ++i){		
		double h = dynamic_cast<FOnePeak*>(m_fp[i])->height();
		if (h> m_maxHeight) m_maxHeight=h;
	}

	m_globalOpt.clear();
	m_numGOpt = 0;
	for (int i = 0; i < m_numPeak; ++i){
		double h = dynamic_cast<FOnePeak*>(m_fp[i])->height();
		if (h == m_maxHeight) {
			m_globalOpt.appendOptima(Solution<CodeVReal>(CodeVReal(move(mapFromPeak(m_fp[i]->getGOpt()[0].data().m_x, m_fp[i]->gidx(), m_fp[i]->idx())), m_fp[i]->getGOpt()[0].data().m_obj)),false);
			++m_numGOpt;
		}
	}


}

void FFreePeak_D_OnePeak::setShapeSet(vector<int> *sset) {
	m_shapeSet.clear();
	if (sset) {
		m_shapeSet = *sset;
	}
	else {//default shapes
		for (int i = 1; i <= 8; ++i) m_shapeSet.push_back(i);
	}
}

void  FFreePeak_D_OnePeak::changeShape() {
	for (size_t i = 0; i < m_changingPeak.size(); ++i) {
		FOnePeak*p = dynamic_cast<FOnePeak*>(m_fp[m_changingPeak[i]]);
		int shape = m_shapeSet[Global::msp_global->getRandInt(0, m_shapeSet.size(),Program_Problem) % m_shapeSet.size()];
		p->setShape(shape);
	}
}