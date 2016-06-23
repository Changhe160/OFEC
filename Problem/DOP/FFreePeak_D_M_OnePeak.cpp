#include "FFreePeak_D_M_OnePeak.h"
#include "../FunctionOpt/FOnePeak.h"
#ifdef OFEC_DEMON
#include "../../../ui/Buffer/Scene.h"
extern unique_ptr<Scene> msp_buffer;
extern bool g_algTermination;
#endif

FreePeak_D_M_OnePeak::FreePeak_D_M_OnePeak(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), v[param_numObj]), \
FreePeak(v), m_radius(m_numBox), m_shiftSeverity(v[param_shiftLength]), m_type(TYPE((int)v[param_changeType])), m_centerPeak(v[param_interTest1]), \
m_onCirclePeak(v[param_interTest2]), m_changeInterval(v[param_changeFre]), m_heightSeverity(v[param_jumpHeight]), m_ipsr(1), m_isPR(m_numBox), m_fcb(m_numBox,true){
	if (m_numDim < 2) throw myException("minimum number of dimensions is two @FreePeak_D_M_OnePeak(ParamMap &v)");
	m_peaksPerBox=m_numObj;
	m_initNumDim = m_numDim;
	m_initNumObj = m_numObj;
	initialize();
	m_tag.insert(DOP);
}

void FreePeak_D_M_OnePeak::createPeaks_(){
	m_numPeak = 0;
	m_fp.clear();
	//m_objIdx.clear();
	m_peakIdx.clear();
	string name;
	if (m_name == "FUN_FreePeak_D_M_OnePeak"){
		name = "FUN_OnePeak";
		//name = "FUN_R_OnePeak";	
	}
	else if (m_name == "FUN_FreePeak_D_M_C_OnePeak"){
		name = "FUN_C_OnePeak";
	}
	for (int t = 0, bidx = 0; t < m_tree.size(); ++t){
		for (int i = 0; i < m_numBox; ++i,++bidx){
			(*m_box)[bidx].peak.resize(m_peaksPerBox);

			double h = mc_maxHeight;
			if (m_type == CT_Jump){
				h = h * Global::msp_global->mp_uniformPro->Next();
				m_radius[bidx] = Global::g_arg[param_radius];
			}
			else{
				m_radius[bidx] = Global::msp_global->mp_uniformPro->Next();

			}

			for (int j = 0; j < m_peaksPerBox; ++j){
				FOnePeak *p = nullptr;
				if (bidx == 0){
					if (j == 0){ // set a center peak					
						p = new FOnePeak(Global::msm_pro[name], m_numDim, name, false, m_centerPeak, h);
						if (m_centerPeak == FOnePeak::SH7 || m_centerPeak == FOnePeak::SH8 || m_centerPeak == FOnePeak::SH9 || m_centerPeak == FOnePeak::SH10) p->setRadiusStepFun(m_radius[i]);
					}
					else{//set surrounding peaks
						MyVector v(m_numDim);
						p = new FOnePeak(Global::msm_pro[name], m_numDim, name, false, m_onCirclePeak);//
						v.randOnRadi(m_radius[i] * p->getNearestDis(), Program_Problem);
						p->setLocation(v.data());
					}
				}
				else{//copy the first region
					vector<double> &loc = dynamic_cast<FOnePeak*>((*m_box)[0].peak[j].get())->location();
					if (j == 0)		p = new FOnePeak(Global::msm_pro[name], m_numDim, name, loc, m_centerPeak, h);
					else p = new FOnePeak(Global::msm_pro[name], m_numDim, name, loc, m_onCirclePeak);
					if (j == 0 && (m_centerPeak == FOnePeak::SH7 || m_centerPeak == FOnePeak::SH8 || m_centerPeak == FOnePeak::SH9 || m_centerPeak == FOnePeak::SH10)) p->setRadiusStepFun(m_radius[i]);
				}
				p->gidx() = i;
				p->idx() = j;
				p->tidx() = t;
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
	update_psr();
	
}

void FreePeak_D_M_OnePeak::update_psr(){
	double min = DBL_MAX, max = DBL_MIN;
	for (int i = 0; i < m_numBox; ++i){
		if (m_type == CT_Jump){
			if (dynamic_cast<FOnePeak*>((*m_box)[i].peak[0].get())->height()>max){
				max = dynamic_cast<FOnePeak*>((*m_box)[i].peak[0].get())->height();
				m_ipsr[0] = i;
				m_isPR[i] = true;
			}
			else m_isPR[i] = false;
		}
		else{
			if (m_radius[i]<min){
				min = m_radius[i];
				m_ipsr[0] = i;
				m_isPR[i] = true;
			}
			else m_isPR[i] = false;

		}
	}
	
}
void FreePeak_D_M_OnePeak::parameterSetting(Problem * rP){

}
ReturnFlag FreePeak_D_M_OnePeak::evaluate_(VirtualEncoding &ss, bool rFlag, ProgramMode mode, bool flag_){
	CodeVReal &s = dynamic_cast< CodeVReal&>(ss);
	int tidx = getTreeRoot(s.m_x);
	int gidx = m_tree[tidx]->get_regionIdx(s.m_x);
	int bidx = boxIdx(tidx, gidx);

	for (size_t i = 0; i < (*m_box)[bidx].peak.size(); ++i){
		CodeVReal x(mapToPeak(s.m_x, gidx, i,tidx), 1);
		(*m_box)[bidx].peak[i]->evaluate_(x, rFlag, mode, false);
		ss.m_obj[i] = x.m_obj[0];
		if (m_tl&&!m_fcb[bidx]){
			pair<int, int> pidx = make_pair(bidx, i);
			double dis = getDistance(s, m_peak[m_peakIdx[pidx]], DIS_EUCLIDEAN);
			if (s.getObjDistance(m_peak[m_peakIdx[pidx]].m_obj) <= m_accuracy&&dis<0.01*(*m_box)[bidx].boxSize) m_fcb[bidx] = true;
		}
	}
	if (rFlag)		m_evals++;

	bool flag;
	#ifdef OFEC_CONSOLE
		if (Global::msp_global->mp_algorithm != nullptr)	flag = !Global::msp_global->mp_algorithm->ifTerminating();
		else flag = true;
	#endif
	#ifdef OFEC_DEMON
		flag = true;
	#endif
	if (rFlag&&m_evals%m_changeInterval == 0 && flag){
		change();
		updatePeak();
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
void FreePeak_D_M_OnePeak::change(){
	++m_counter;
	switch (m_type)
	{
	case FreePeak_D_M_OnePeak::CT_Web:
		web();
		break;
	case FreePeak_D_M_OnePeak::CT_Rotation:
		rotate();
		break;
	case FreePeak_D_M_OnePeak::CT_Jump:
		jump();
		break;
	case FreePeak_D_M_OnePeak::CT_Basin:
		varyBasin();
		break;
	case FreePeak_D_M_OnePeak::CT_NumObj:
		varyNumObj();
		break;
	case FreePeak_D_M_OnePeak::CT_NumDim:
		varyNumDim();
		break;
	}
	update_psr();

	for (auto i = 0; i < m_fcb.size(); ++i){
		if (m_tl)	m_fcb[i] = false;
		else m_fcb[i] = true;
	}
}
void FreePeak_D_M_OnePeak::jump(){
	for (int i = 0; i < m_numBox; ++i){
		if (!m_fcb[i]) continue;
		FOnePeak*p = dynamic_cast<FOnePeak*>((*m_box)[i].peak[0].get());
		double h = p->height();
		double delta = 7*m_heightSeverity*Global::msp_global->mp_normalPro->Next();
		if (h + delta < 0) h -= delta;
		else h += delta;
		if (h>mc_maxHeight) h = mc_maxHeight*Global::msp_global->mp_uniformPro->Next();
		p->changeHeight(h);
	}
}
void FreePeak_D_M_OnePeak::rotate(){

	vector<int> d(m_numDim);
	Matrix I(m_numDim, m_numDim), rotM(m_numDim, m_numDim);
	double angle = -OFEC_PI + 2 * OFEC_PI*m_shiftSeverity;
	for (int i = 0; i<m_numBox; i++){
		if (!m_fcb[i]) continue;
		Global::msp_global->initializeRandomArray(d, m_numDim, Program_Problem);
		I.identity();		
		for (int j = 0; j + 1<m_numDim; j += 2){			
			I.setRotationAngle(d[j], d[j + 1], angle);
			if (j == 0) rotM = I;
			else		rotM *= I;
		}
		Matrix m(m_numDim, 1);

		for (int j = 0; j < m_peaksPerBox; ++j){
			FOnePeak*p = dynamic_cast<FOnePeak*>((*m_box)[i].peak[j].get());
			vector<double> loc =p->location();
			m.setDataRow(loc.data(), m_numDim);
			m *= rotM;
			copy(m[0].begin(), m[0].end(), loc.begin());
			p->shiftLocation(loc);
		}
	}
	
}

void FreePeak_D_M_OnePeak::web(){
	for (int i = 0; i<m_numBox; i++){
		if (!m_fcb[i]) continue;
		double delta = m_shiftSeverity*(Global::msp_global->mp_uniformPro->Next()-0.5);
		double radi = m_radius[i];
		if (radi + delta < 0 || radi + delta>1){
			radi -= delta;
		}
		else{
			radi += delta;
		}
		m_radius[i] = radi;
		FOnePeak* pc = dynamic_cast<FOnePeak*>((*m_box)[i].peak[0].get());
		MyVector center(pc->location());
		for (int j = 1; j < m_peaksPerBox; ++j){
			FOnePeak*p = dynamic_cast<FOnePeak*>((*m_box)[i].peak[j].get());
			MyVector loc = p->location();
			loc -= center;
			loc.normalize();
			loc *= m_radius[i]*pc->getNearestDis();
			p->shiftLocation(loc.data());
		}
	}
}
void FreePeak_D_M_OnePeak::varyBasin(){
	BI b;
	m_basin.clear();
	set<int> ts;
	for (int i = 0; i < m_numBox; ++i){
		int tidx, gidx;
		treeIdx(i, tidx, gidx);
		ts.insert(tidx);
		m_tree[tidx]->get_leafParent(gidx, b.idx, b.cutdim, b.low, b.high);
		vector<BI>::iterator it = m_basin.begin();
		for (; it != m_basin.end(); ++it){
			if (it->idx == b.idx) break; // already in m_basin
		}
		if (it == m_basin.end()) m_basin.push_back(b);
	}

	for (auto it = m_basin.begin(); it != m_basin.end(); ++it){
		double mag = (it->high - it->low)*m_shiftSeverity/2;
		double delta = Global::msp_global->getRandFloat(-mag, mag,Program_Problem);
		int tidx, gidx;
		treeIdx(it->idx, tidx, gidx);

		double val = m_divisionPoint[tidx][gidx][it->cutdim];
		if (val + delta>it->high || val + delta < it->low){
			val -= delta;
		}
		else{
			val += delta;
		}
		m_divisionPoint[tidx][gidx][it->cutdim] = val;
	}
	for (auto i:ts)	m_tree[i]->buildIndex();
	computeBoxSize();
}
void FreePeak_D_M_OnePeak::varyNumDim(){
	
	int newdims;
	int period = 25, offset = m_initNumDim;
	int minDim = 2;
	if ((m_counter + offset - minDim) / period % 2 == 0) newdims = ((m_counter + offset - minDim) % period);
	else newdims = period - ((m_counter + offset - minDim) % 25);

	newdims += 2; //minimum number of dimensions is two
	Global::g_arg[param_numDim] = newdims;
	FreePeak_D_M_OnePeak *newpro = new FreePeak_D_M_OnePeak(Global::g_arg);
	Global::g_arg[param_numDim] = m_initNumDim;

	if (newdims < m_numDim){
		vector<int> cd(newdims);
		vector<int> sd(m_numDim);
		Global::msp_global->initializeRandomArray(sd, Program_Problem);
		copy(sd.begin(), sd.begin() + newdims, cd.begin());
		newpro->copyChanges(this,&cd);
	}else	newpro->copyChanges(this);

	FreePeak::resizeDim(newdims);
	*this = *newpro;
	delete newpro;
}
void FreePeak_D_M_OnePeak::varyNumObj(){

	int newobjs;
	int period = 25, offset = m_initNumObj;
	int min = 2;
	if ((m_counter + offset - min) / period % 2 == 0) newobjs = ((m_counter + offset - min) % period);
	else newobjs = period - ((m_counter + offset - min) % 25);

	newobjs += 2; //minimum number of dimensions is two
	Global::g_arg[param_numObj] = newobjs;
	FreePeak_D_M_OnePeak *newpro = new FreePeak_D_M_OnePeak(Global::g_arg);
	Global::g_arg[param_numObj] = m_initNumObj;

	if (newobjs < m_numObj){
		vector<int> co(newobjs);
		vector<int> so(m_numObj);
		Global::msp_global->initializeRandomArray(so, Program_Problem);
		copy(so.begin(), so.begin() + newobjs, co.begin());
		newpro->copyChanges(this, nullptr,&co);
	}
	else	newpro->copyChanges(this);

	newpro->copyChanges(this);
	FreePeak::resizeObj(newobjs);
	*this = *newpro;
	delete newpro;

}

void FreePeak_D_M_OnePeak::copyChanges(const Problem * op, const vector<int> *cd, const vector<int> *co){

	FreePeak::copyChanges(op,cd,co);
	const FreePeak_D_M_OnePeak *pro = dynamic_cast<const FreePeak_D_M_OnePeak*>(op);
	m_radius = pro->m_radius;
	m_counter = pro->m_counter;    // the current number of changes
	m_initNumDim = pro->m_initNumDim;
	m_initNumObj=pro->m_initNumObj;

}

void FreePeak_D_M_OnePeak::varyNoBoxes(){
// TODO
}