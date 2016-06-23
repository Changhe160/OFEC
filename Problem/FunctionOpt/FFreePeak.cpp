#include "FFreePeak.h"

FreePeak::FreePeak(ParamMap &v) : ContinuousProblem((v[param_proId]), (v[param_numDim]), (v[param_proName]), v[param_numObj]), m_box(new vector<Box>(v[param_numBox])), m_numBox(v[param_numBox])
{
	if (v.find(param_peaksPerBox) != v.end()) 
		m_peaksPerBox = v[param_peaksPerBox];
	else m_peaksPerBox =1;

	if (v.find(param_divisionMode) != v.end()) m_dm = v[param_divisionMode];

	vector<string> s = { "","Random","Even","UserDefined" };
	m_proPar << "Number of boxes: " << m_numBox << "; Division: " << s[m_dm] << ";";
}

void FreePeak::initialize(){
	set<ProTag> tag;
	tag.insert(CONT);
	if (m_numObj > 1) tag.insert(MOP);
	else tag.insert(SOP);
	setProTag(tag);
	setAccuracy(1.e-5);
	setSearchRange(-100, 100);
	
	vector<vector<pair<double, double>>> box;

	if (Global::g_arg.find(param_dataFile2) != Global::g_arg.end()){
		stringstream ss;
		ss << Global::g_arg[param_workingDir] << Global::g_arg[param_dataDirectory1] << Global::g_arg[param_dataFile2];
		ifstream in(ss.str().c_str());
		if (!in) throw myException("file does not exist @ FreePeak::initialize()");
		int tnum,nnum,tint;
		string line,temp;
		getline(in, line);
		istringstream iss(line);
		iss >>tnum >>nnum;
		m_numBox = nnum;
		m_box.reset(new vector<Box>(m_numBox));
		
		for (int t = 0; t < tnum; ++t){
			getline(in, line);
			istringstream iss(line);
			int node;
			iss >> temp >> tint >> node;
			m_divisionPoint.push_back(vector<vector<double>>(node, vector<double>(m_numDim)));
			getline(in, line);
			iss.clear();
			iss.str(line);
			
			box.push_back(move(vector<pair<double, double>>(m_numDim)));
			char c;
			for (int d = 0; d < m_numDim; ++d){
				iss >>c>> box.back()[d].first >>c>> box.back()[d].second>>c;
			}
			for (int i = 0; i < node; ++i){
				iss.clear();
				getline(in, line);
				iss.str(line);
				for (int d = 0; d < m_numDim; ++d)  iss >> m_divisionPoint[t][i][d];
			}
		}
		in.close();
	}
	else{
		box.push_back(move(vector<pair<double, double>>(m_numDim,pair<double,double>(-100,100))));
		m_divisionPoint.push_back(vector<vector<double>>(m_numBox - 1, vector<double>(m_numDim)));
		if (m_dm == DivisionMode::DM_Random){
			for (int i = 0; i < m_numBox - 1; ++i){
				for (auto &j : m_divisionPoint[0][i]){
					j = -100 + 200 * Global::msp_global->mp_uniformPro->Next();
				}
			}
		}
		else if (m_dm == DivisionMode::DM_Even){
			int div = 0;
			vector<pair<double, double>> ran(m_numDim, pair<double, double>(-100.0, 100.0));
			for (int i = 0; i < m_numBox - 1; ++i){
				double r = 1. / (m_numBox - i);
				div = i%m_numDim;
				for (int j = 0; j < m_numDim; ++j){
					if (j == div){
						m_divisionPoint[0][i][j] = ran[j].first + (ran[j].second - ran[j].first)*r;
						ran[j].first = m_divisionPoint[0][i][j];
					}
					else{
						m_divisionPoint[0][i][j] = (ran[j].first + ran[j].second) / 2;
					}
				}
			}
		}
	}

	for (auto i = 0; i < box.size(); ++i){
		m_tree.push_back(unique_ptr<KDTreeSpace::PartitioningKDTree<double>>(new KDTreeSpace::PartitioningKDTree<double>(m_numDim, m_divisionPoint[i])));
		m_tree[i]->setInitBox(box[i]);
		m_tree[i]->buildIndex(m_dm);
	}

	computeBoxSize();
	createPeaks();
}
void FreePeak::setBoxSize(const vector<double>& rat){

	for (int t = 0; t < m_tree.size(); ++t){
		int div = 0;
		m_tree[t]->get_rootBox();
		vector<pair<double, double>> ran(m_tree[t]->get_rootBox().box);
		double ratio = 1;
		for (int j = 0; j < m_numDim; ++j) ratio *= 200/(ran[j].second - ran[j].first) ;
		for (int i = 0; i < m_divisionPoint[t].size(); ++i){
			double r = rat[i] * ratio;
			div = i%m_numDim;
			for (int j = 0; j < m_numDim; ++j){
				if (j == div){
					m_divisionPoint[t][i][j] = ran[j].first + (ran[j].second - ran[j].first)*r;
					ran[j].first = m_divisionPoint[t][i][j];
				}
				else{
					m_divisionPoint[t][i][j] = (ran[j].first + ran[j].second) / 2;
				}
			}
		}

		m_tree[t]->buildIndex(m_dm);
	}

	computeBoxSize();

}
void FreePeak::createPeaks(){
	//a subresion of m_tree->region[i] associates with m_peak[i]
	createPeaks_();
	m_peak.resize(m_numPeak, CodeVReal(m_numDim, m_numObj));
	updatePeak();
	setDisAccuracy(disNearestPeak() / 2);
	setOptType((*m_box)[0].peak[0]->getOptType(),-1); // assume all subproblems have the same opt. type
	updateGOpt();
}


double FreePeak::disNearestPeak(){
	double min =DBL_MAX;
	for (int i = 0; i < m_numPeak; ++i){
		for (int j = 0; j < i; ++j){
			double dis = getDistance(m_peak[i], m_peak[j],DIS_EUCLIDEAN);
			if (dis<min) min = dis;
		}
	}
	return min;
}
void FreePeak::computeBoxSize(){
	int idx = 0;
	for (int t = 0; t < m_tree.size(); ++t){
		double r = 1.,l,u;
		for (int i = 0; i < m_numDim; ++i){
			m_searchRange.getSearchRange(l, u, i);
			r *= (m_tree[t]->get_rootBox().box[i].second - m_tree[t]->get_rootBox().box[i].first) / (u - l);
		}
		for (int i = 0; i < m_tree[t]->region.size(); ++i,++idx){
			(*m_box)[idx].boxRatio = m_tree[t]->region[i].rat*r;
			(*m_box)[idx].boxSize = (*m_box)[idx].boxRatio* m_searchRange.getDomainSize();
		}
	}
}
void FreePeak::updateGOpt(){
	if (m_numObj == 1){
		m_globalOpt.setFlagLocTrue();
		Solution<CodeVReal> opt((*m_box)[0].peak[0]->getGOpt()[0]);
		for (int i = 1; i < m_numPeak; ++i){
			if (compare(opt.data(), m_fp[i]->getGOpt()[0].data()) == Compare_Worse)opt = m_fp[i]->getGOpt()[0];
		}

		m_globalOpt.clear();
		for (int i = 0; i < m_numBox; ++i){
			CompareResultFlag flag =compare( opt.data(),(*m_box)[i].peak[0]->getGOpt()[0].data());
			if (flag == Compare_Worse || flag == Compare_Equal){
				for (int j = 0; j < (*m_box)[i].peak[0]->getGOpt().getNumOpt(); ++j){
					m_globalOpt.appendOptima(Solution<CodeVReal>(move(CodeVReal(move(mapFromPeak((*m_box)[i].peak[0]->getGOpt()[j].data().m_x, (*m_box)[i].peak[0]->gidx(), (*m_box)[i].peak[0]->idx(), (*m_box)[i].peak[0]->tidx())), (*m_box)[i].peak[0]->getGOpt()[j].data().m_obj))), false);
				}
			}
		}
	}
}

void FreePeak::updatePeak(){
	for (int p = 0; p < m_numPeak; ++p){
		CodeVReal peak(m_numDim, m_numObj);
		peak.m_x = move(mapFromPeak(m_fp[p]->getGOpt()[0].data().m_x, m_fp[p]->gidx(), m_fp[p]->idx(), m_fp[p]->tidx()));
		evaluate_(peak, false, Program_Problem, false);
		m_peak[p] = peak;
	}
}
void FreePeak::copyChanges(const Problem * op, const vector<int> *cd, const vector<int> *co){
	ContinuousProblem::copyChanges(op,cd,co);

	const FreePeak *pro = dynamic_cast<const FreePeak*>(op);
	for (int t = 0; t < m_tree.size(); ++t){
		for (int i = 0; i < m_tree[t]->region.size(); ++i){
			if (cd){
				for (int j = 0; j < m_numDim; ++j){
					m_divisionPoint[t][i][j] = pro->m_divisionPoint[t][i][(*cd)[j]];
				}
			}
			else		copy(pro->m_divisionPoint[t][i].begin(), pro->m_divisionPoint[t][i].end(), m_divisionPoint[t][i].begin());
		}
		m_tree[t]->buildIndex();
	}
	
	computeBoxSize();

	for (int i = 0; i < m_numBox; ++i){
		if (co){
			for (int j = 0; j < m_peaksPerBox; ++j){
				(*m_box)[i].peak[j]->copyChanges((*pro->m_box)[i].peak[(*co)[j]].get(), cd);
			}
		}
		else{
			for (int j = 0; j < m_peaksPerBox; ++j){
				(*m_box)[i].peak[j]->copyChanges((*pro->m_box)[i].peak[j].get(), cd);
			}
		}
	}

	updatePeak();

}

FreePeak& FreePeak::operator = (FreePeak& rhs){
	ContinuousProblem::operator=(rhs);

	m_box.reset(rhs.m_box.release());	
	m_numPeak = rhs.m_numPeak;
	m_numBox=rhs.m_numBox;
	m_divisionPoint=rhs.m_divisionPoint;
	for (int t = 0; t < m_tree.size(); ++t){
		m_tree[t].reset(rhs.m_tree[t].release());
	}
	m_peaksPerBox = rhs.m_peaksPerBox;
	m_fp=rhs.m_fp;
	m_objIdx;
	m_peakIdx=rhs.m_peakIdx;
	m_peak=rhs.m_peak;
	return *this;
}

void  FreePeak::resizeDim(int num){
	ContinuousProblem::resizeDim(num);
	for (auto&i : m_divisionPoint){
		for (auto &j : i) j.resize(num);
	} 
	for (auto &i : m_peak) i.m_x.resize(num);


}
void  FreePeak::resizeObj(int num){
	ContinuousProblem::resizeObj(num);
	for (auto &i : m_peak) i.m_obj.resize(num);
}

void FreePeak::initializeSolution(VirtualEncoding &s, const int idx, const int maxIdx){
	CodeVReal&result = dynamic_cast< CodeVReal&>(s);
	if (m_initialBox4Sol.size() == 0){
		ContinuousProblem::initializeSolution(s, idx, maxIdx);
	}
	else{
		// randomly initialize solotions in allowed boxes
		int idx = Global::msp_global->getRandInt(0, m_initialBox4Sol.size(), Program_Problem);
		double l, u;
		int bidx = m_initialBox4Sol[idx];
		int tidx, gidx;
		treeIdx(bidx, tidx, gidx);
		for (int i = 0; i<m_numDim; i++){
			u = m_tree[tidx]->region[gidx].box[i].second;
			l = m_tree[tidx]->region[gidx].box[i].first;
			result[i] = l + (u - l)*Global::msp_global->mp_uniformAlg->Next();
		}
	}
}

void FreePeak::generateDivision(const vector<vector<pair<double, double>>> &box, vector<int> &node, int numNode, int numDim, int numTree, const string &path, bool evendiv, vector<double> *br) {
	stringstream ss;
	static int fileno = 1;
	ss <<"div"<< fileno << "_dim" << numDim << "_tree" << numTree << "_node" << numNode << ".div";

	string file = path +ss.str();
	ofstream out(file);

	Global::g_arg[param_dataFile2] = ss.str();
	out << node.size() << " " << numNode << endl;

	if (br != 0|| evendiv) {
		int bidx = 0;
		for (int t = 0; t < numTree; ++t) {
			out << "Tree " << t + 1 << " " << node[t] << endl;
			for (int d = 0; d < numDim; ++d) out << "[ " << box[t][d].first << " : " << box[t][d].second << " ] ";
			out << endl;
			vector<pair<double, double>> ran(box[t]);

			for (int n = 0; n < node[t]; ++n,++bidx) {
				double x,r;
				int div;
				if (br) r = (*br)[bidx];
				else 	r = 1. / (node[t] - n);

				div = n%numDim;
				for (int j = 0; j < numDim; ++j) {
					if (j == div) {
						x = ran[j].first + (ran[j].second - ran[j].first)*r;
						ran[j].first = x;
					}
					else {
						x = (ran[j].first + ran[j].second) / 2;
					}
					out << x << " ";
				}
				out << endl;
			}
		}
	}else {

		default_random_engine gen;
		uniform_real_distribution<double> udis(0.0, 1.0);
		for (int t = 0; t < numTree; ++t) {
			out << "Tree " << t + 1 << " " << node[t] << endl;
			for (int d = 0; d < numDim; ++d) out << "[ " << box[t][d].first << " : " << box[t][d].second << " ] ";
			out << endl;
			for (int n = 0; n < node[t]; ++n) {
				double x;
				for (int d = 0; d < numDim; ++d) {
					x = box[t][d].first + (box[t][d].second - box[t][d].first)*udis(gen);
					out << x << " ";
				}
				out << endl;
			}
		}
	}	
	out << "#END";
	out.close();
	fileno++;
}

void FreePeak::generateLocationSingleObj(int numDim, int numPeak, const string & path, bool random) {

	static int fileno = 1;

	stringstream ss;
	ss << "location" << fileno<< "_dim" << numDim << "_peak" << numPeak << ".loc";
	string file = path + ss.str();
	ofstream out(file);
	Global::g_arg[param_dataFile3] = ss.str();

	out << "#PEAK \t DIM" << endl;
	out << numPeak << " " << numDim << endl;
	out << "#LOCATION" << endl;

	if (random) {
		default_random_engine gen;
		uniform_real_distribution<double> udis(0.0, 1.0);
		double x;
		for (int jj = 0; jj < numPeak; ++jj) {
			for (int d = 0; d < numDim; ++d) {
				x = -100 + 200 * udis(gen);
				out << x << " ";
			}
			out << endl;
		}
	}
	else {
		for (int jj = 0; jj < numPeak; ++jj) {
			for (int d = 0; d < numDim; ++d) {
				out << " 0 ";
			}
			out << endl;
		}
	}

	out << "#END";
	out.close();
	fileno++;
}

void FreePeak::generateLocationMultiObj(int numDim, int numBox, int numObj, float minRadius, const string & path, int mode,float rhoweb) {
	static int fileno = 1;
	vector<vector<double>> peak(numObj, vector<double>(numDim, 0));
	default_random_engine gen;
	uniform_real_distribution<double> udis(0.0, 1.0);

	for (int oi = 1; oi < numObj; ++oi) {
		MyVector v(numDim);
		v.randOnRadi(minRadius * 100, udis, gen);
		peak[oi] = v.data();
	}

	stringstream ss;
	ss << "location" << fileno<< "_dim" << numDim << "_obj" << numObj << "_box" << numBox << ".loc";
	string file = path + ss.str();
	ofstream out(file);
	Global::g_arg[param_dataFile3] = ss.str();

	out << "#BOX \t OBJ \t DIM" << endl;
	out << numBox << " " << numObj << " " << numDim << endl;
	out << "#LOCATION \t RADIUS" << endl;
	if (mode == 1) { //Jump case	
		for (int pb = 0; pb < numBox; ++pb) {
			for (int p = 0; p < numObj; ++p) {
				for (int d = 0; d < numDim; ++d) {
					out << peak[p][d] << " ";
				}
				if (p == 0) out << minRadius;
				else out << 0;
				out << endl;
			}
		}	
	}
	else if(mode==2){//Web case	
		vector<double> radius(numObj, 0);
		radius[0] = minRadius;

		for (int pb = 0; pb < numBox; ++pb) {
			if (pb != 0) {
				radius[0] = minRadius + rhoweb*udis(gen)*(1 - minRadius);
				for (int oi = 1; oi < numObj; ++oi) {
					MyVector v(numDim);
					v.randOnRadi(100 * radius[0], udis, gen);
					peak[oi] = v.data();
				}
			}
			for (int p = 0; p < numObj; ++p) {
				for (int d = 0; d < numDim; ++d) {
					out << peak[p][d] << " ";
				}
				out << radius[p];
				out << endl;
			}
		}
	}
	else if (mode == 3) {// Countable
		for (int pb = 0; pb < numBox; ++pb) {
			MyVector v(numDim);
			v.randomize(udis,gen,-100,100);			
			for (int p = 0; p < numObj; ++p) {
				for (int d = 0; d < numDim; ++d) {
					out << v[d] << " ";
				}
				out << 0<<endl;
			}
		}
	}
		
	out << "#END";
	out.close();
	fileno++;
}

void FreePeak::generateConfGOP(int numPeak, const string & path,int gShape, int gTrans, const string & gPosition, int lShape, int lTrans, double lMaxH, bool lRandH,bool trap, double tMinH, double tMaxH) {
	vector<Peak> conf(numPeak);	
	default_random_engine gen;
	uniform_real_distribution<double> udis(0.0, 1.0);
	vector<int> gbox;
	vector<int> shape = { 1,2,3,4,5,7,8,9,10,11,12,13 };
	for (int p = 0; p < numPeak; ++p) {
		int curuse;
		if (trap) { // set trap at the first box
			if(p==0){
			conf[p].basin = "0";
			conf[p].height = tMaxH;
			conf[p].minHeigh = tMinH;
			conf[p].shape = 6;
			conf[p].transf.push_back(0); //without transformation
			}
			else if (p == 1) {//global optimum
				conf[p].basin = gPosition;
				conf[p].height = 100;
				conf[p].minHeigh = 0;
				conf[p].shape = gShape;
				conf[p].transf.push_back(gTrans);
				gbox.push_back(p);
			}
			else {	//local optimimum			
				conf[p].basin = "random";
				if (lRandH) conf[p].height = lMaxH * udis(gen);
				else	conf[p].height = lMaxH;
				conf[p].minHeigh = 0;
				if (lShape == -1) {
					conf[p].shape = shape[udis(gen)*shape.size()];
				}else	conf[p].shape = lShape;
				conf[p].transf.push_back(lTrans);
			}
		}else if(p==0){//global optimum
			conf[p].basin = "0";
			conf[p].height = 100;
			conf[p].minHeigh = 0;
			conf[p].shape = gShape;
			conf[p].transf.push_back(gTrans);
			gbox.push_back(p);
		}
		else {	//local optimimum			
			conf[p].basin = "random";
			if(lRandH) conf[p].height = lMaxH * udis(gen);
			else	conf[p].height = lMaxH;
			conf[p].minHeigh = 0;
			if (lShape == -1) {
				conf[p].shape = shape[udis(gen)*shape.size()];
			}else	conf[p].shape = lShape;
			conf[p].transf.push_back(lTrans);
		}		
	}


	generateConf(conf, gbox,path, 1);

}

void  FreePeak::generateConf(const vector<Peak> &conf, const vector<int> &gbox, const string & path, int obj) {
	static int fileno = 1;

	stringstream ss;
	ss << "config" << fileno << "_box" << conf.size() / obj << "_obj" << obj << ".conf";
	string file = path + ss.str();
	ofstream out(file);
	Global::g_arg[param_dataFile1] = ss.str();

	out << "#BOX \t OBJ " << endl;
	out << conf.size() / obj << " " << obj << endl;
	out << "#BASIN OF PARATO/GLOBAL OPT." << endl;
	for (auto i = 0; i < gbox.size()-1; ++i) out << gbox[i] << ":";
	out << gbox.back() << endl;

	out << "#BASIN \t SHAPE \t HEIGHT \t MINHEIGHT \t TRANSFORMATION" << endl;
	for (auto i = 0; i < conf.size(); ++i) {
		out << conf[i].basin << " " << conf[i].shape << " " << conf[i].height << " " << conf[i].minHeigh << " ";
		for (auto j = 0; j < conf[i].transf.size(); ++j) {
			if (j < conf[i].transf.size() - 1) 	out << conf[i].transf[j] << ":";
			else out << conf[i].transf[j];
		}
		out << endl;
	}
	out << "#END" << endl;
	out.close();
	fileno++;
}

void FreePeak::generateConfMMOP(int numPeak, const string & path, int gShape, int gTrans, const string & gPosition, int gNum, int lShape, int lTrans, double lMaxH, bool lRandH, bool trap, int tNum, double tMinH, double tMaxH) {
	vector<Peak> conf(numPeak);
	default_random_engine gen;
	uniform_real_distribution<double> urd(0.0, 1.0);
	vector<int> gshape = { 1, 2, 6,11, 12, 13 };
	vector<int> lshape = { 1,2,3,4,5,7,8,9,10,11,12,13 };
	vector<int> gbox;
	for (int p = 0; p < numPeak; ++p) {
		if (p <tNum) {
			conf[p].basin = "random";
			conf[p].height = tMaxH;
			conf[p].minHeigh = tMinH;
			conf[p].shape = 6;
			conf[p].transf.push_back(0); //without transformation
		}
		else if (p <tNum+ gNum) {//global optimum
			conf[p].basin = "random";
			conf[p].height = 100;
			conf[p].minHeigh = 0;
			if (gShape == -1)	conf[p].shape = gshape[urd(gen)*gshape.size()];
			else	conf[p].shape = gShape;
			conf[p].transf.push_back(gTrans);
			gbox.push_back(p);
		}
		else {	//local
			conf[p].basin ="random";
			if (lRandH) conf[p].height = lMaxH * urd(gen);
			else	conf[p].height = lMaxH;
			conf[p].minHeigh = 0;
			
			if (lShape == -1) {
				conf[p].shape = lshape[urd(gen)*lshape.size()];
			}
			else conf[p].shape = lShape;

			conf[p].transf.push_back(0);
		}
	}
	generateConf(conf, gbox,path, 1);
}

void FreePeak::generateConfMOP(int numBox, const string & path, int o1stShape, int numObj, int mode, float rhoHeight,int numParetoRegion) {

	vector<Peak> conf(numBox*numObj);
	vector<int> bidx(numBox);
	for (int p = 0; p < numBox; ++p)bidx[p] = (p);
	default_random_engine gen;
	uniform_real_distribution<double> udis(0.0, 1.0);
	vector<int> gbox;
	if (mode = 1 || mode == 2) { //Jump or Web
		gbox.push_back(0);
		for (int bi = 0; bi < numBox; ++bi) {
			int ridx;
			if (bi == 0) ridx = 0;
			else ridx = udis(gen)*bidx.size();
			for (int oidx = 0; oidx < numObj; oidx++) {
				int p = bi*numObj + oidx;
				if (p / numObj == 0) {//pareto region
					conf[p].basin = "0";
					conf[p].height = 100;
				}
				else {	//non-pareto region
					conf[p].basin = to_string(bidx[ridx]);
					if (mode == 1 && p%numObj == 0) {
						conf[p].height = 100 * (1 - udis(gen)*rhoHeight);
					}
					else {
						conf[p].height = 100;
					}
				}
				conf[p].minHeigh = 0;
				if (oidx == 0)	conf[p].shape = o1stShape;
				else conf[p].shape = 1;
				conf[p].transf.push_back(0);
			}

			auto it = bidx.begin();
			while (*it != bidx[ridx]) ++it;
			bidx.erase(it);
		}
	}
	else if(mode==3){ //Countable
		vector<vector<double>> obj(numParetoRegion, vector<double>(numObj));
		for (int i = 0; i < numParetoRegion; ++i) {
			for (int j = 0; j < numObj; ++j) {
				if (j == 0) {
					obj[i][j] = 100 - i *20. / numParetoRegion;
				}
				else if (j == 1) {
					obj[i][j] = 80 + (i + 1) *20. / numParetoRegion;
				}
				else {
					obj[i][j] = 100* udis(gen);
				}
			}
		}
		vector<int> idx(numBox);
		for (int i = 0; i < numBox; ++i) idx[i] = i;
		for (int i = 0; i < numParetoRegion; ++i){
			int j = udis(gen)*idx.size();
			gbox.push_back(idx[j]);		
			idx.erase(idx.begin() + j);
		}
		std::sort(gbox.begin(), gbox.end());
		vector<int> shape = { 1, 2, 3,4,5,6,7,8,9,10,11, 12, 13 };
		for (int bi = 0,j=0; bi < numBox; ++bi) {		
			for (int oidx = 0; oidx < numObj; oidx++) {
				int p = bi*numObj + oidx;
				conf[p].basin = to_string(bi);

				if (bi == gbox[j]) {
					conf[p].height = obj[j][oidx];
				}
				
				conf[p].minHeigh = 0;
				
				conf[p].shape = shape[shape.size()*udis(gen)];
				conf[p].transf.push_back(0);
			}
			if (bi == gbox[j]) ++j;			
		}
	}
	
	generateConf(conf, gbox, path, numObj);

}

void  FreePeak::setup(const string & problem) {
	//default settings

	int numDim = 5, numObj = 1, numBox = 2;
	Global::g_arg[param_numDim] = numDim;
	Global::g_arg[param_numBox] = numBox;
	Global::g_arg[param_numObj] = numObj;
	vector<int> node;
	typedef  pair<double, double> RANGE;
	vector<vector<RANGE>> box(1, vector<RANGE>(numDim, RANGE(-100, 100)));
	node.push_back(numBox - 1);

	string path;
	
	if (problem == "FUN_FreePeak_OnePeak") {
		int gShape = 2;/*[1,13]*/ int gTrans = 0; /*[0,5] */ const string gPosition = "random";/* {  "first","smallest", "largest", "random" }*/
		int lShape = -1; /*[-1,1-13]*/int lTrans = 0; double lMaxH = 99; bool lRandH = true;
		bool trap = true; double tMinH = 90; double tMaxH = 95;

		/*path = (string) Global::g_arg[param_workingDir];
		path += "Problem/data/GOP/";
		Global::g_arg[param_dataDirectory1]= "Problem/data/GOP/";
		generateDivision(box, node, numBox, numDim, 1, path);
		generateLocationSingleObj(numDim, numBox, path, false);
		generateConfGOP(numBox, path, gShape, gTrans, gPosition, lShape, lTrans, lMaxH, lRandH, trap, tMinH, tMaxH);
		*/
		path = (string)Global::g_arg[param_workingDir];
		path += "Problem/data/GOP/";
		Global::g_arg[param_dataDirectory1] = "Problem/data/GOP/";
		vector<double> br(numBox);
		default_random_engine gen;
		uniform_real_distribution<double> udis(0.0, 1.0);
		double sum = 0;
		for (auto&i : br) {
			i = udis(gen);
			sum += i;
		}

		for (auto&i : br) i /= sum;
		generateDivision(box, node, numBox, numDim, 1, path, false, &br);
		generateLocationSingleObj(numDim, numBox, path, false);
		generateConfGOP(numBox, path, gShape, gTrans, gPosition, lShape, lTrans, lMaxH, lRandH, trap, tMinH, tMaxH);
		/*path = (string) Global::g_arg[param_workingDir];
		path += "Problem/data/MMOP/";
		Global::g_arg[param_dataDirectory1] = "Problem/data/MMOP/";
		int gNum = numBox*0.4,tNum=numBox*0.4;
		gShape = -1;
		generateDivision(box, node, numBox, numDim, 1, path);
		generateLocationSingleObj(numDim, numBox, path, false);
		generateConfMMOP(numBox, path, gShape, gTrans, gPosition, gNum,lShape, lTrans, lMaxH, lRandH, trap,tNum, tMinH, tMaxH);*/
	}
	else if (problem == "FUN_FreePeak_D_OnePeak") {
	
	}
	else if (problem == "FUN_FreePeak_M_OnePeak") {
		numObj = 2;
		Global::g_arg[param_numObj] = numObj;
		int type = 1;
		Global::g_arg[param_case] = type;

		int o1stShape = 1;/*[1,8]*/ float rhoHeight = 0.1; /*[0,1]*/
		int numParetoRegion = 1;
		if (type == 3) numParetoRegion = 0.5*numBox;

		Global::g_arg[param_numParetoRegion] = numParetoRegion;

		path = (string)Global::g_arg[param_workingDir];
		path += "Problem/data/MOP/";
		Global::g_arg[param_dataDirectory1] = "Problem/data/MOP/";
		generateDivision(box, node, numBox, numDim, 1, path);
		generateLocationMultiObj(numDim, numBox, numObj, 0.1, path, type);
		generateConfMOP(numBox, path, o1stShape, numObj, type, rhoHeight, numParetoRegion);
	}
	else if (problem == "FUN_FreePeak_D_M_OnePeak"){

	}

}