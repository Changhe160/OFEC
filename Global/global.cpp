#include "global.h"

#ifdef OFEC_CONSOLE
thread_local unique_ptr<Global> Global::msp_global=nullptr; 
#endif

#ifdef OFEC_DEMON
unique_ptr<Global>  Global::msp_global;
#endif

STRING2ID Global::msm_alg,Global::msm_pro, Global::msm_term;
map<ALG2PRO,unsigned> Global::msm_alg4pro;
classFactory  Global::ms_classFactory;
map<string,Param> Global::msm_param;
int Global::ms_curProId,Global::ms_curAlgId;
ParamMap Global::g_arg;

int &g_curProID=Global::ms_curProId;

#ifdef OFEC_DEMON
extern bool g_algTermination;
#endif

Global::Global(int runId,double seedPro, double seedAlg):mp_problem(nullptr),mp_algorithm(nullptr),m_totalNumIndis(0),m_runId(runId){
	double seed=0;
	if(seedPro<0||seedPro>1) seed=0.5;
	else seed=seedPro;
	 mp_cauchyPro=move(unique_ptr<Cauchy>(new Cauchy(seed)));
	 mp_normalPro=move(unique_ptr<Normal>(new Normal(seed)));
	 mp_uniformPro=move(unique_ptr<Uniform>(new Uniform(seed)));
	 mp_levyPro=move(unique_ptr<Levy>(new Levy(1.4,seed)));

	 if(seedAlg<0||seedAlg>1) seed=0.5;
	 else seed=seedAlg;

	 mp_cauchyAlg=move(unique_ptr<Cauchy>(new Cauchy(seed)));
	 mp_normalAlg=move(unique_ptr<Normal>(new Normal(seed)));
	 mp_uniformAlg=move(unique_ptr<Uniform>(new Uniform(seed)));
	 mp_levyAlg=move(unique_ptr<Levy>(new Levy(1.4,seed)));

}

Global::~Global(){	 
	 mp_cauchyPro.reset();
	 mp_normalPro.reset();
	 mp_uniformPro.reset();
	 mp_levyPro.reset();
	 mp_cauchyAlg.reset();
	 mp_normalAlg.reset();
	 mp_uniformAlg.reset();
	 mp_levyAlg.reset();
	 mp_problem.reset();
	 mp_algorithm.reset();
}

void Global::setSeedAlg(double seed){

	 mp_cauchyAlg.reset(new Cauchy(seed));
	 mp_normalAlg.reset(new Normal(seed));
	 mp_uniformAlg.reset(new Uniform(seed));
	 mp_levyAlg.reset(new Levy(1.4,seed));

}

int Global::getRandInt(const int min, const int max,ProgramMode mode){
	if (min == max||min==max-1) return min;
	if(mode==Program_Algorithm)	return static_cast<int>(min+(max-min)*mp_uniformAlg->Next());
	else	return static_cast<int>(min+(max-min)*mp_uniformPro->Next());
}
double Global::getRandFloat(const double min, const double max,ProgramMode mode){
	if (min == max) return min;
	if(mode==Program_Algorithm)	return min+(max-min)*mp_uniformAlg->Next();
	else		return min+(max-min)*mp_uniformPro->Next();
}


bool Global::registerAlg4Pro(ALG2PRO  alg2pro){
	return registerItem(msm_alg4pro,alg2pro);
}
void Global::registerParamter(){
	msm_param["param_ND"]=param_numDim;
	msm_param["param_NP"]=param_numPeak;
	msm_param["param_PN"]=param_proName;
	msm_param["param_AN"]=param_algName;
	msm_param["param_ME"]=param_maxEvals;
	msm_param["param_SL"]=param_shiftLength;
	msm_param["param_CT"]=param_changeType;
	msm_param["param_CR"]=param_changeRatio;
	msm_param["param_RID"]=param_runId;
	msm_param["param_AID"]=param_algId;
	msm_param["param_PID"]=param_proId;
	msm_param["param_FNDC"]=param_flagNumDimChange;
	msm_param["param_FNPC"]=param_flagNumPeakChange;
	msm_param["param_PNCM"]=param_peakNumChangeMode; 
	msm_param["param_FN"]=param_flagNoise;
	msm_param["param_FTL"]=param_flagTimeLinkage;
	msm_param["param_CDBGFID"]=param_comDBGFunID;
	msm_param["param_NS"]=param_noiseSeverity;
	msm_param["param_TLS"]=param_timelinkageSeverity;
	msm_param["param_PS"]=param_popSize;
	msm_param["param_ECF"]=param_evalCountFlag;
	msm_param["param_SI"]=param_stepIndi;
	msm_param["param_TT"]=param_trainingTime;
	msm_param["param_SPS"]=param_subPopSize;
	msm_param["param_OLD"]=param_overlapDgre;
	msm_param["param_CTH"]=param_convThreshold;
	msm_param["param_TW"]=param_timeWindow;
	msm_param["param_SF"]=param_sampleFre;
	msm_param["param_gOptFlag"]=param_gOptFlag;
	msm_param["param_WD"]=param_workingDir;
	msm_param["param_CF"]=param_changeFre;
	msm_param["param_CoF"]=param_convFactor;
	msm_param["param_NR"]=param_numRun;
	msm_param["param_NT"]=param_numTask;
	msm_param["param_ER"]=param_exlRadius;
	msm_param["param_SVM"]=param_solutionValidationMode;
	msm_param["param_PIM"]=param_populationInitialMethod;
	msm_param["param_HR"]=param_hibernatingRadius;
	msm_param["param_MNPS"]=param_minNumPopSize;
	msm_param["param_NC"]=param_numChange;
	msm_param["param_NO"]=param_numObj;
	msm_param["param_C"]=param_case;
	msm_param["param_RBP"]=param_resource4BestPop;
	msm_param["param_DF1"]=param_dataFile1;
	msm_param["param_XP"]=param_xoverProbability;
	msm_param["param_MP"]=param_mutProbability;
	msm_param["param_PT"] = param_proTag;
	msm_param["param_NGO"] = param_numGOpt;
	msm_param["param_NF"] = param_noiseFlag;
	msm_param["param_PF"] = param_predicFlag;
	msm_param["param_CT2"] = param_changeType2;
	msm_param["param_PPB"] = param_peaksPerBox;
	msm_param["param_IT1"]=param_interTest1;
	msm_param["param_IT2"]=param_interTest2;
	msm_param["param_IT3"]=param_interTest3;
	msm_param["param_IT4"]=param_interTest4;
	msm_param["param_NB"] = param_numBox;
	msm_param["param_HCM"] = param_heightConfigMode;
	msm_param["param_PC"] = param_peakCenter;
	msm_param["param_NPR"] = param_numParetoRegion;
	msm_param["param_VR"] = param_validRadius;
	msm_param["param_A"] = param_attraction;
	msm_param["param_R"] = param_radius;
	msm_param["param_JH"] = param_jumpHeight; 
	msm_param["param_VR"] = param_variableRelation; 
	msm_param["param_PkS"] = param_peakShape; 
	msm_param["param_DM"] = param_divisionMode; 
	msm_param["param_POS"] = param_peakOffset;
	msm_param["param_FI"] = param_flagIrregular;
	msm_param["param_FA"] = param_flagAsymmetry;
	msm_param["param_DF2"] = param_dataFile2;
	msm_param["param_DF3"] = param_dataFile3;
	msm_param["param_FR"] = param_flagRotation;
	msm_param["param_DD1"] = param_dataDirectory1;
	msm_param["param_MI"] = param_maxIter;
	msm_param["param_MSI"] = param_maxSucIter;
	msm_param["param_E"] = param_epsilon;
	msm_param["param_MSDE"] = param_mutationSchemeDE;
	msm_param["param_USPL"] = param_updateSchemeProbabilityLearning;
}

string gGetProblemName(const int id){
	STRING2ID::iterator it;
	for(it=Global::msm_pro.begin();it!=Global::msm_pro.end();it++){ 
		if(it->second==id) return it->first;
	}
	return "";
}
string gGetAlgorithmName(const int id){
	STRING2ID::iterator it;
	for(it=Global::msm_alg.begin();it!=Global::msm_alg.end();it++){ 
		if(it->second==id) return it->first;
	}
	return "";
}


char *gStrtok_r(char *s, const char *delim, char **save_ptr){   

    if (s == 0) s = *save_ptr;   
   
    /* Scan leading delimiters.  */   
    s += strspn(s, delim);   
    if (*s == '\0') return 0;   
   
    /* Find the end of the token.  */   
    char *token = s;   
    s = strpbrk(token, delim);   
    if (s == 0)   
        /* This token finishes the string.  */   
        *save_ptr = strchr(token, '\0');   
    else{   
        /* Terminate the token and make *SAVE_PTR point past it.  */   
        *s = '\0';   
        *save_ptr = s + 1;   
    }   
   
    return token;   
}  

bool gIsDynamicProlem(){
	string pro = gGetProblemName(Global::ms_curProId);

	if (pro == "FUN_FreePeak_D_OnePeak" || pro == "FUN_FreePeak_D_R_OnePeak" || pro == "FUN_FreePeak_D_C_OnePeak" || pro == "FUN_FreePeak_D_M_OnePeak"\
		|| pro == "DYN_CONT_CompositionDBG" || pro == "DYN_CONT_RotationDBG" || pro == "DYN_CONT_MovingPeak") return true;
	else return false;
}

void gGetStringValue(const string &str, vector<TypeVar>& val){
	//val=1:2:3
	string s(str);
	string letter;
	for (auto i = 'a'; i <= 'z'; i++) {
		letter += i;
		letter += i - 32;
	}
	size_t pos1 = s.find("="), pos2;
	s.erase(0, pos1 + 1);
	val.clear();
	while (1){
		pos2 = s.find(":");
		string value;
		if (pos2 != string::npos){
			value = s.substr(0, pos2);
		}
		else{
			value = s.substr(0, s.size());
		}

		if (value.compare("true") == 0){
			val.push_back(true);
		}
		else if (value.compare("false") == 0){
			val.push_back(false);
		}
		else if (value.find_first_of(letter) != string::npos){
			if (value.size() == 1)
				val.push_back(value[0]);
			else
				val.push_back(value);
		}
		else if (string::npos != value.find('.')){
			val.push_back(atof(value.c_str()));
		}
		else{
			val.push_back(atoi(value.c_str()));
		}
		if (pos2 != string::npos)	s.erase(0, pos2 + 1);
		else break;
	}
}

std::istream& gSafeGetline(std::istream& is, std::string& t){
	t.clear();

	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}