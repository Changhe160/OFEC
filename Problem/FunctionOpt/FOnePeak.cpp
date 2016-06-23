#include "FOnePeak.h"

FOnePeak::FOnePeak(ParamMap &v) :Problem((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1), \
BenchmarkFunction((v[param_proId]), (v[param_numDim]), (v[param_proName]), 1), m_location(m_numDim), m_transfLoc(m_numDim), m_vr(v[param_variableRelation]), \
m_dominoWeight(v[param_numDim]),m_irr(v[param_flagIrregular]),m_asy(v[param_flagAsymmetry]){
	setSearchRange(-100, 100);
	m_rotationFlag = v[param_flagRotation];
	m_shape = v[param_peakShape];
	initialize();
	if (v[param_peakCenter] == 1){
		setLocationAtCenter();
	}
	else if (v[param_peakCenter] == 2){
		setLocationRandom();
	}
	
	//difficulty();
}
FOnePeak::FOnePeak(const int rId, const int rDimNumber, string& rName, bool random, int shape, double h, vector<int> *transform, int label ) :Problem(rId, rDimNumber, rName, 1), \
BenchmarkFunction(rId, rDimNumber, rName, 1), m_location(m_numDim), m_shape(shape), m_height(h), m_dominoWeight(rDimNumber), m_label(label), m_transfLoc(m_numDim){
	setSearchRange(-100, 100);
	setFlag(transform);	
	initialize();
	if (!random)	setLocationAtCenter();
	else setLocationRandom();
}
FOnePeak::FOnePeak(const int rId, const int rDimNumber, string& rName, const vector<double>&loc, int shape, double h, vector<int> *transform, int label) :Problem(rId, rDimNumber, rName, 1), \
BenchmarkFunction(rId, rDimNumber, rName, 1), m_location(m_numDim), m_shape(shape), m_height(h), m_dominoWeight(rDimNumber), m_label(label),m_transfLoc(m_numDim){
	setSearchRange(-100, 100);
	setFlag(transform);
	initialize();
	setLocation(loc);
}
void FOnePeak::setFlag(vector<int> *trans){
	if (trans){
		m_vr = (*trans)[0];
		m_rotationFlag = (*trans)[1];	
		m_irr = (*trans)[2];
		m_asy = (*trans)[3];		
	}
	else{
		m_rotationFlag = false;
		m_irr = false;
		m_asy = false;
		m_vr = Separable;
	}
}
void FOnePeak::initialize(){
	if (m_rotationFlag){
		setConditionNumber(2);
		loadRotation();
	}
	setOptType(MAX_OPT);
	setDisAccuracy(0.2);
	setAccuracy(1.e-6);
	m_globalOpt.setFlagLocTrue();
	m_globalOpt[0].data().m_obj[0] = m_height;
	m_originalGlobalOpt = m_globalOpt;
	if (m_shape == SH8 || m_shape == SH7 || m_shape == SH9 || m_shape == SH10)	initilizeStepFun();
	if (m_vr == FOnePeak::VariableRelat::Domino){
		for (auto&i : m_dominoWeight) i = pow(10, Global::msp_global->mp_normalPro->Next());
	}
	if (m_vr == FOnePeak::VariableRelat::PartSeparable) {
		m_noGroup = Global::msp_global->getRandInt(2, m_numDim, Program_Problem);
		m_group.resize(m_noGroup);
		vector<int> set(m_numDim-1),randIdx(m_numDim);
		for (int i = 0; i < m_numDim-1; ++i) {
			set[i] = i;
		}
		Global::msp_global->initializeRandomArray(randIdx, m_numDim);
		vector<int> pick(m_noGroup);
		for (int i = 0; i < m_noGroup-1; ++i) {
			pick[i] = Global::msp_global->randPick(set,Program_Problem);
		}
		pick.back() = m_numDim - 1;
		std::sort(pick.begin(), pick.end());

		for (int i = 0,j=0; i < m_noGroup; ++i) {
			while (j <= pick[i]) {
				m_group[i].push_back(randIdx[j++]);
			}
		}
	}
	setObjSet();
}
void FOnePeak::initilizeStepFun(bool flag){	
	if (Global::g_arg.find(param_radius) != Global::g_arg.end()){
		setRadiusStepFun(Global::g_arg[param_radius],flag);
	}
	else{
		setRadiusStepFun(0.5,flag);
	}
}
void FOnePeak::setStepFun(){
	m_radius = m_ratio*m_nearestDis;
	if (m_shape == SH9)	m_foldHeight = -m_height / (sqrt(m_radius + 1));
	else if (m_shape == SH8) m_foldHeight = m_height*(-1 - m_theta*(m_k-1)) / (sqrt(m_radius + 1));
	else m_foldHeight = -m_height; 
}
void FOnePeak::setRadiusStepFun(double r, bool flag){
	m_ratio = r;	
	if (flag) setStepFun();
}

void FOnePeak::setLocationAtCenter(){ //default setup
	for (int i = 0; i < m_numDim; ++i){
		m_location[i] = (m_searchRange[i].m_lower + m_searchRange[i].m_upper) / 2;
	}
	m_transfLoc = m_location;
	transform(m_transfLoc.data());
	m_globalOpt[0].data().m_x = m_location;
	m_originalGlobalOpt[0].data().m_x = m_location;
	computeDisToBoarder();
	computeMinHeight();
}
void FOnePeak::setLocationRandom(){
	for (int i = 0; i < m_numDim; ++i){
		m_location[i] = m_searchRange[i].m_lower + (m_searchRange[i].m_upper - m_searchRange[i].m_lower)*Global::msp_global->mp_uniformPro->Next();
	}
	m_transfLoc = m_location;
	transform(m_transfLoc.data());
	m_globalOpt[0].data().m_x = m_location;
	m_originalGlobalOpt[0].data().m_x = m_location;
	computeDisToBoarder();
	computeMinHeight();
}
void FOnePeak::setLocation(const vector<double> &x){
	m_location = x;
	m_transfLoc = m_location;
	transform(m_transfLoc.data());
	m_globalOpt[0].data().m_x = m_location;
	m_originalGlobalOpt[0].data().m_x = m_location;
	computeDisToBoarder();
	computeMinHeight();
}
void FOnePeak::evaluate__(double const *x_, vector<double>& obj){
	double dummy = getTransDistance(x_);

	switch (m_shape)
	{
	case SH1: obj[0] = m_height - dummy;		break;	//linear
	case SH2: obj[0] = m_height*exp(-dummy); break;	//convex
	case SH3: obj[0] = m_height - sqrt(m_height)*sqrt(dummy); break;//covex
	case SH4: obj[0] = m_height / (dummy + 1); break;	//convex
	case SH5: obj[0] = m_height - dummy*dummy / m_height; break;//concave
	case SH6: obj[0] = m_height - exp(2 * sqrt(dummy / sqrt(m_numDim))) + 1; break;//concave	
	case SH7: //concave+linear
		if (dummy<m_radius)		obj[0] = m_height * cos(OFEC_PI*dummy / m_radius);
		else obj[0] = m_foldHeight - (dummy - m_radius);
		break;	
	case SH8:{
		if (dummy < m_radius){
			int domain = floor(m_k*dummy / m_radius);
			obj[0] = m_height*(cos(m_k * OFEC_PI*(dummy - domain*m_radius / m_k) / m_radius) - m_theta*domain) / (sqrt(dummy + 1));
		}
		else{
			obj[0] = m_foldHeight - (dummy - m_radius);
		}
	}	
		break;
	case SH9: //concave+convex+linear
		if (dummy <= m_radius){
			obj[0] = m_height*cos((2*m_k-1) * OFEC_PI*dummy / m_radius) / (sqrt(dummy+1));
		}
		else{
			obj[0] = m_foldHeight - (dummy - m_radius);
		}		
		break; 	
	case SH10:	//flat, nightmare when radius=0
		if (dummy <= m_radius)		obj[0] = m_height;
		else obj[0] = m_foldHeight ;
		break;	
	case SH11:	//netural+needle
		if (dummy <= 10)		obj[0] = m_height*exp(-dummy);
		else obj[0] = m_height*exp(-10);
		break;	
	case SH12:	//ripple
		obj[0] = m_height*cos((2 * m_k - 1) * OFEC_PI*dummy / m_furestDis) / sqrt(log(dummy + 1)+1);
		break;
	case SH13:	//netural+needle+ripple
		if (dummy <= 10)		obj[0] = m_height*exp(-dummy);
		else if (dummy>10 && dummy<=50) obj[0] = m_height*exp(-10);
		else obj[0] = m_height*sin((2 * m_k - 1) * OFEC_PI*(dummy - 50) / (m_furestDis - 50)) / sqrt(log(dummy - 49) + 1) + m_height*exp(-10);
		break;
	}

	if (m_standardize) obj[0] = m_stdMinHeight + (m_stdMaxHeight - m_stdMinHeight)*(obj[0] - m_minHeight) / (m_height - m_minHeight);
}
void FOnePeak::setHeight(const double h){ 
	m_height = h; 
	m_globalOpt[0].data().m_obj[0] = m_height;
	m_originalGlobalOpt[0].data().m_obj[0] = m_height;
	if (m_shape == SH8 || m_shape == SH7 || m_shape == SH9 || m_shape == SH10)	setStepFun();
	computeMinHeight();
}
bool FOnePeak::isValid(const VirtualEncoding  &s){
	if (m_name == "FUN_C_OnePeak"){
		//a region within a sphere with radius r is valid 
		const CodeVReal  &x = dynamic_cast<const CodeVReal  &>(s);
		double dis = 0;
		for (int i = 0; i < m_numDim; ++i){
			dis += (m_location[i] - x[i])*(m_location[i] - x[i]);
		}
		dis = sqrt(dis);
		if (dis > m_searchRange.getDomainSize()*m_validRadius) return false;
		else return true;		
	}
	else{
		return ContinuousProblem::isValid(s);
	}
}
void FOnePeak::setValidRadius(double r){
	m_validRadius = r;
}
void FOnePeak::shiftLocation(const vector<double> &x){
	if (m_maxMemorysize > 0){
		m_preLoc.push_front(m_location);
		if (m_preLoc.size() > m_maxMemorysize) m_preLoc.pop_back();
	}
	m_location = x;
	SolutionValidation mode = VALIDATION_REMAP;
	validate_(m_location, &mode);
	m_globalOpt[0].data().m_x = m_location;
	m_originalGlobalOpt[0].data().m_x = m_location;
	computeDisToBoarder();
	computeMinHeight();
}
void FOnePeak::shiftLocation(){
	if (m_maxMemorysize > 0){
		m_preLoc.push_front(m_location);
		if (m_preLoc.size() > m_maxMemorysize) m_preLoc.pop_back();
	}
	
	MyVector v(m_numDim);
	v.randOnRadi(m_shiftSeverity*m_searchRange.getDomainSize(), Program_Problem);
	for (int i = 0; i < m_numDim; ++i){
		m_location[i] += v[i];
	}	
	SolutionValidation mode = VALIDATION_REMAP;
	validate_(m_location, &mode);
	m_globalOpt[0].data().m_x = m_location;
	m_originalGlobalOpt[0].data().m_x = m_location;
	computeDisToBoarder();
	computeMinHeight();
}

void FOnePeak::changeHeight(const double h){
	if (m_maxMemorysize > 0){
		m_preHeight.push_front(m_height);
		if (m_preHeight.size() > m_maxMemorysize) m_preHeight.pop_back();
	}
	setHeight(h);
	updateRobust();
	m_quality = m_height / m_maxHeight;
}
void FOnePeak::setMemorySize(const int val){
	m_maxMemorysize = val;
}
void FOnePeak::setMaxHeight(const double val){
	m_maxHeight = val;
	m_quality = m_height / m_maxHeight;
}
void  FOnePeak::updateRobust(){
	
	if (m_preHeight.size() > 0){
		double delta = 0;
		delta += fabs(m_height - m_preHeight[0]);
		
		for (auto i = m_preHeight.begin(); i != m_preHeight.end() - 1; ++i){
			delta += fabs(*i - *(i + 1));
		}	
		delta /= m_preHeight.size();
		m_robustness = delta;
	}
	else{
		m_robustness = 0;
	}
	
}

void FOnePeak::computeMinHeight(){
	double dummy = m_furestDis;
	switch (m_shape)
	{
	case SH1: m_minHeight = m_height - dummy; break;
	case SH2: m_minHeight = m_height*exp(-dummy); break;
	case SH3: m_minHeight = m_height - sqrt(m_height)*sqrt(dummy); break;
	case SH4: m_minHeight = m_height / (dummy + 1); break;
	case SH5: m_minHeight = m_height - dummy*dummy / m_height; break;
	case SH6: m_minHeight = m_height - exp(2*sqrt(dummy/sqrt(m_numDim))) + 1; break;
	case SH7: case SH8:
		m_minHeight = m_foldHeight - (dummy - m_radius);
		break;
	case SH9:{
		double minInR = -m_height / (sqrt(m_radius / (2 * m_k - 1) + 1));
		dummy = m_foldHeight - (dummy - m_radius);
		m_minHeight = minInR<dummy ? minInR : dummy;
	}
		   break; 
	case SH10:
		m_minHeight = m_foldHeight;
		break;
	case SH11:
		m_minHeight = m_height*exp(-10);
		break;
	case SH12:
		m_minHeight = -m_height / sqrt(log(m_furestDis/(2*m_k-1) + 1) + 1);
		break;
	case SH13:
		m_minHeight = -m_height / sqrt(log(1.5*(m_furestDis - 50) / (2 * m_k - 1) +1) + 1) + m_height*exp(-10);
		break;
	}
}

void FOnePeak::computeDisToBoarder(){
	double dummy = 0;
	m_nearestDis = DBL_MAX;
	for (int i = 0; i < m_numDim; ++i){
		double l, u;
		m_searchRange.getSearchRange(l, u, i);
		if (u - m_location[i] > m_location[i] - l){
			dummy += (u - m_location[i])*(u - m_location[i]);
			if (m_nearestDis > m_location[i] - l) m_nearestDis = m_location[i] - l;
		}
		else{
			dummy += (m_location[i] - l)*(m_location[i] - l);
			if (m_nearestDis > u - m_location[i])m_nearestDis = u - m_location[i];
		}
	}
	m_furestDis = sqrt(dummy);
	if (m_shape == SH8 || m_shape == SH7 || m_shape == SH9 || m_shape == SH10)	setStepFun();
}

void FOnePeak::difficulty(double rh){
	
	vector<double> sample(10000);
	static bool flag=false;
	static vector<CodeVReal> points(10000, CodeVReal(m_numDim, m_numObj));
	for (int i = 0; i < 10000; ++i){		
		if (!flag){
			for (auto &j : points[i].m_x){
				j = -100 + Global::msp_global->mp_uniformPro->Next() * 200;
			}
		}		
		evaluate__(points[i].m_x.data(), points[i].m_obj);
		sample[i] = points[i].m_obj[0];
	}
	flag = true;
	double max, min;
	max = *max_element(sample.begin(), sample.end());
	min = *min_element(sample.begin(), sample.end());
	for_each(sample.begin(), sample.end(), [&](double &n){n = rh*(n - min) / (max - min); });
	
	m_mean = accumulate(sample.begin(), sample.end(), 0.0);
	m_mean /= sample.size();

	vector<double> diff(sample.size());
	std::transform(sample.begin(), sample.end(), diff.begin(), std::bind2nd(std::minus<double>(), m_mean));
	m_variance = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / (diff.size()-1);
	
}
void FOnePeak::setBasinRatio(double r){
	m_basinRatio = r;
}
void FOnePeak::copyChanges(const Problem * pro, const vector<int> *cd, const vector<int> *co){

	const FOnePeak *op = dynamic_cast<const FOnePeak *>(pro);

	if (cd){
		for (int j = 0; j < m_numDim; ++j){
			m_location[j] = op->m_location[(*cd)[j]];
			m_transfLoc[j] = op->m_transfLoc[(*cd)[j]];
		}
	}
	else{
		copy(op->m_location.begin(), op->m_location.end(), m_location.begin());
		copy(op->m_transfLoc.begin(), op->m_transfLoc.end(), m_transfLoc.begin());
	}
	
	computeDisToBoarder();
	computeMinHeight();

	m_preHeight=op->m_preHeight;
	m_preLoc=op->m_preLoc;		// take care of # of dimensions
	m_height = op->m_height;
	m_shape = op->m_shape;
	m_standardize = op->m_standardize;
	m_robustness=op->m_robustness;
	m_maxMemorysize = op->m_maxMemorysize;
	m_shiftSeverity = op->m_shiftSeverity;
	m_heightSeverity=op->m_heightSeverity;
	m_quality=op->m_quality;
	m_maxHeight = op->m_maxHeight;
	m_validRadius = op->m_validRadius;
	m_radius=op->m_radius;
	m_foldHeight=op->m_foldHeight;
	m_ratio =op->m_ratio;
	m_basinRatio = op->m_basinRatio;
	m_dominoWeight = op->m_dominoWeight;
	m_k = op->m_k;
	m_theta = op->m_theta;

	m_globalOpt[0].data().m_x = m_location;
	m_globalOpt[0].data().m_obj[0] = m_height;
	m_originalGlobalOpt = m_globalOpt;

	m_stdMaxHeight = op->m_stdMaxHeight;
	m_stdMinHeight = op->m_stdMinHeight;
	m_label = op->m_label;
}

void FOnePeak::irregularize(double *x){
	double c1, c2, a;
	for (int i = 0; i < m_numDim; ++i){
		if (x[i]>0){
			c1 = 10;	c2 = 7.9;
		}else {
			c1 = 5.5;	c2 = 3.1;
		}
		if (x[i] != 0){
			a = log(fabs(x[i]));
		}
		else a = 0;
		x[i] = gSign(x[i])*exp(a + 0.049*(sin(c1*x[i] + sin(c2*x[i]))));
	}
}
void FOnePeak::asyemmetricalize(double *x){
	double belta = 0.1;
	if (m_numDim == 1) return;
	
	for (int i = 0; i < m_numDim; ++i){
		if (x[i]>0){
			x[i] = pow(x[i], 1. + belta*(i / (m_numDim - 1))*log(x[i]+1));
		}
	}
		
}

bool FOnePeak::loadRotation(){
	string s;
	stringstream ss;
	ss << m_numDim << "Dim.txt";
	s = ss.str();
	s.insert(0, m_name + "_Rot_");

	s.insert(0, "Problem/FunctionOpt/Data_PFs/GOP/");
	s.insert(0, Global::g_arg[param_workingDir]);
	ifstream in;
	in.open(s.data());
	if (in.fail()){		
		mp_rotationMatrix->randomize(Program_Problem);
		mp_rotationMatrix->generateRotationMatrix(m_conditionNumber, Program_Problem);
		
		ofstream out(s.c_str());
		mp_rotationMatrix->Print(out);
		out.close();
	}
	else{
		mp_rotationMatrix->Read_Data(in);
	}
	in.close();

	return true;
}
void FOnePeak::transform(double * x){

	double *x_ = new double[m_numDim];	
	for (int i = 0; i < m_numDim; ++i) x_[i] = x[i] - m_location[i];
	
	if (m_irr)		irregularize(x_);
	if (m_asy)		asyemmetricalize(x_);
	if (m_rotationFlag){
		for (int i = 0; i<m_numDim; i++) {
			x[i] = 0;

			for (int j = 0; j < m_numDim; j++) {
				x[i] += (*mp_rotationMatrix)[j][i] * x_[j];
			}
		}
	}
	else copy(x_, x_ + m_numDim, x);

	delete[] x_;
}

vector<double> FOnePeak::getPoint(double oval, const MyVector& dir, int ith){

	vector<double> p;
	MyVector v(dir);
	double dummy = getDistance2Center(oval,ith);
	if (dummy >= 0){
		p.resize(m_numDim);
		v.normalize();
		v *= dummy;
		copy(v.data().begin(), v.data().end(), p.begin());
	}
	return move(p);
}

double FOnePeak::getDistance2Center(double oval,int ith){
	bool flag = m_irr || m_asy || m_rotationFlag || m_vr != FOnePeak::VariableRelat::Separable;
	if (flag) return -1;
	if (m_standardize){
		oval = (oval - m_stdMinHeight)*(m_height - m_minHeight) / (m_stdMaxHeight - m_stdMinHeight) + m_minHeight;
	}  

	double dummy = -1;
	switch (m_shape){
	case SH1: dummy = m_height - oval; 		break;
	case SH2: dummy = -log(oval / m_height); break;
	case SH3: dummy = pow((m_height - oval) / sqrt(m_height), 2.); break;
	case SH4: dummy = (m_height / oval) - 1; break;
	case SH5:  dummy = sqrt((m_height - oval)*m_height); break;
	case SH6: dummy = pow(log(m_height - oval + 1) / 2, 2.)*sqrt(m_numDim); break;
	case SH7:
		if (oval >= m_foldHeight){
			dummy = m_radius*acos(oval / m_height) / OFEC_PI;
		}
		else{
			dummy = m_foldHeight - oval + m_radius;
		}
		break;
	case SH8:
		if (oval > m_foldHeight){
			double d1, d2, val, d3;
			d1 = (ith - 1)*m_radius / m_k, d2 = ith*m_radius / m_k;
			do{
				d3 = (d1 + d2) / 2;
				val = m_height*(cos(m_k * OFEC_PI*(d3 - floor(m_k*d3 / m_radius)) / m_radius) - m_theta*floor(m_k*d3 / m_radius)) / (sqrt(d3 + 1));
				if (val > oval){
					d1 = d3;
				}
				else{
					d2 = d3;
				}
			} while (fabs(val - oval)>1.e-10);
			dummy = d3;
		}
		else{
			dummy = m_foldHeight - oval + m_radius;
		}
		break;
	case SH9:{		
		double d1 = 0, d2 = m_radius / (2 * m_k - 1), val, d3;
		d1 = (ith - 1)*m_radius / (2 * m_k - 1), d2 = ith*m_radius / (2 * m_k - 1);
		do{
			d3 = (d1 + d2) / 2;
			val = m_height*cos((2 * m_k - 1) * OFEC_PI*d3 / m_radius) / (sqrt(d3 + 1));
			if (val > oval){
				d1 = d3;
			}
			else{
				d2 = d3;
			}
		} while (fabs(val - oval)>1.e-10);
		dummy = d3;		
	}
		break;
	default:
		break;
	}
	
	return dummy;
}

void FOnePeak::getContiSegment(const vector<double> pos, vector<pair<vector<double>, vector<double>>> &seg){
	seg.clear();
	CodeVReal x(m_numDim, m_numObj);
	switch (m_shape){
	case SH1: 	case SH2:	case SH3:	case SH4: 	case SH5: 	case SH6: case SH7:
		seg.push_back(make_pair(m_location, pos));
		break;
	case SH8:{
		double dummy = getTransDistance(pos.data());
		double step = m_radius / m_k;
		MyVector vn(m_numDim);
		for (int d = 0; d < m_numDim; ++d) vn[d] = pos[d] - m_location[d];
		vn.normalize();
		for (int i = 1; i <=m_k; ++i){
			MyVector v1(vn), v2(vn);
			v1 = vn*(step*(i - 1) + 1.0e-10);
			if (dummy > step*i + 1.0e-10){
				v2 = vn*(step*i - 1.0e-10);
				seg.push_back(make_pair(v1.data(), v2.data()));
			}
			else{								
				seg.push_back(make_pair(v1.data(), pos));
				break;
			}	
		}
		
	}
		break;
	case SH9:{	
		double dummy = getTransDistance(pos.data());
		int i;
		MyVector vn(m_numDim);
		for (int d = 0; d < m_numDim; ++d) vn[d] = pos[d] - m_location[d];
		vn.normalize();
		vector<MyVector> vv;
		
		for (i = 1; i <= m_k; ++i){			
			if (i == 1){
				vv.push_back(MyVector(m_location));
			}
			else{
				MyVector v(vn);
				v *= m_radius*(i - 1) * 2 / (2 * m_k - 1);
				vv.push_back(v);
			}
			if (i < m_k){
				double d = m_radius*i * 2 / (2 * m_k - 1);
				MyVector v(vn);
				v *= d;
				x.m_x = v.data();
				evaluate_(x, false, Program_Problem, false);
				vv.push_back(move(MyVector(getPoint(x.m_obj[0], vn,i*2-1))));
			}
			else{
				double d = m_radius*(i * 2-1) / (2 * m_k - 1);
				MyVector v(vn);
				v *= d;
				vv.push_back(v);
			}
		}

		for (i = 1; i <= m_k; ++i){
			double dis = m_radius*(2 * i - 1) / (2 * m_k - 1);
			if (dummy >= dis) {
				seg.push_back(move(make_pair(vv[2 * i - 2].data(), vv[2 * i - 1].data())));
			}
			else break;
		}
		if (i <= m_k){
			int is = (2 * m_k - 1)*dummy / m_radius;
			if (is % 2 == 0)	seg.push_back(move(make_pair(vv[2 * i - 2].data(), pos)));
			else{
				MyVector mid(m_numDim);
				for (int d = 0; d < m_numDim; ++d) mid[d] = (vv[is + 1][d] + vv[is][d]) / 2;
				double dis = getTransDistance(mid.data().data());

				if (dis<=dummy)	seg.push_back(move(make_pair(pos, vv[is + 1].data())));
				else{
					MyVector x1(pos), x2(vv[is + 1]),x3,x4;
					do{
						mid = x1.getPointBetween(x2, 0.5);
						x.m_x = mid.data();
						evaluate_(x, false, Program_Problem, false);
						x3.data()=getPoint(x.m_obj[0], mid, is);
						x4=mid.getPointBetween(x3, 0.5);
						dis = getTransDistance(x4.data().data());
						if (dis <= dummy){
							x1 = mid;
						}
						else{
							x2 = mid;
						}
					} while (fabs(dis-dummy)>1.0e-10);
					seg.back().second = x3.data();
					seg.push_back(move(make_pair(pos, mid.data())));
				}
			}
		}
		else seg.back().second = pos;
	}
		break;
	default:
		break;
	}
}

void FOnePeak::getContiSegment(const vector<double> p1, const vector<double> p2, vector<pair<vector<double>, vector<double>>> &seg){
	seg.clear();
	CodeVReal x(m_numDim, m_numObj);
	switch (m_shape){
	case SH1: 	case SH2:	case SH3:	case SH4: 	case SH5: 	case SH6: case SH7:{
		seg.push_back(make_pair(p1, p2));
		}	
		break;
	case SH8:{
		MyVector mid(m_numDim), v1(m_numDim), v2(m_numDim), vn1(p1), vn2(p2);
		for (int d = 0; d < m_numDim; ++d) mid[d] = (p1[d] + p2[d])/2;
		double step = m_radius / m_k;
		double dummy = getTransDistance(mid.data().data());
		deque<MyVector> lseg;

		vn1 -= mid;	vn2 -= mid;
		vn1.normalize();	vn2.normalize();
		
		for (int i = 1; i < m_k; i++){
			if (dummy <= step*i){
				v1 = vn1;	v2 = vn2;
				double dis2mid = sqrt(pow(step*i, 2) - dummy*dummy);
				v1 *= dis2mid;	v2 *= dis2mid;
				v1 += mid;		v2 += mid;
				lseg.push_front(v1);
				lseg.push_back(v2);
			}
		}

		lseg.push_front(move(MyVector(p1)));
		lseg.push_back(move(MyVector(p2)));

		int midx = lseg.size() / 2-1;
		for (int i = 0; i < lseg.size()-1; ++i){
			MyVector dir1(lseg[i]), dir2(lseg[i + 1]);
			dir1 -= m_location;
			dir2 -= m_location;
			if (i == midx){
				dir1 *= 1 - 1.e-10;
				dir2 *= 1 - 1.e-10;
			}
			else if(i<midx){
				dir1 *= 1 - 1.e-10;
				dir2 *= 1 + 1.e-10;
			}
			else{
				dir1 *= 1 + 1.e-10;
				dir2 *= 1 - 1.e-10;
			}
			dir1 += m_location;
			dir2 += m_location;
			seg.push_back(make_pair(dir1.data(), dir2.data()));
		}	
	}
		break;
	
	default:
		break;
	}
}

void FOnePeak::setShape(int shape) {
	if (m_shape != shape) {
		m_shape = shape;		
		if (m_shape == SH8 || m_shape == SH7 || m_shape == SH9 || m_shape == SH10)	initilizeStepFun(true);
		computeMinHeight();
	}
}