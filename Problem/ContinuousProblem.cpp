#include "ContinuousProblem.h"
#include "../Utility/myVector.h"
#include "../Global/global.h"

bool ContinuousProblem::isValid(const VirtualEncoding &ss){
	const CodeVReal&s=dynamic_cast<const CodeVReal&>(ss);
	for(int i=0;i<m_numDim;i++){
		if(m_searchRange[i].m_flag){
			if((s[i])<m_searchRange[i].m_lower||(s[i])>m_searchRange[i].m_upper) return false;
		}
	}
	return true;
	
}
void ContinuousProblem::validate_(vector<double> &s, SolutionValidation *mode){
	double l, u, x;
	SolutionValidation val;
	if (mode) val = *mode;
	else val = m_validationMode;
	switch (val){
	case VALIDATION_IGNORE:
		break;
	case VALIDATION_REINITIALIZE:
		for (int i = 0; i<m_numDim; i++){
			m_searchRange.getSearchRange(l, u, i);
			s[i] = l + (u - l)*Global::msp_global->mp_uniformAlg->Next();
		}
		break;
	case VALIDATION_REMAP:
		for (int i = 0; i<m_numDim; i++){
			m_searchRange.getSearchRange(l, u, i);
			x = s[i];
			if (x<l){
				s[i] = l + (u - l)*(l - x) / (u - x);
			}
			else if (x>u){
				s[i] = l + (u - l)*(u - l) / (x - l);
			}
		}
		break;
	case VALIDATION_SETTOBOUND:
		for (int i = 0; i<m_numDim; i++){
			m_searchRange.getSearchRange(l, u, i);
			x = s[i];
			if (x<l){
				s[i] = l;
			}
			else if (x>u){
				s[i] = u;
			}
		}
		break;
	}
}
void ContinuousProblem::validate(VirtualEncoding &ss,SolutionValidation *mode){
	validate_(dynamic_cast< CodeVReal&>(ss).m_x, mode);
}
void ContinuousProblem::initializePartSolution(VirtualEncoding &result,int begin,int end){
	if(begin<0||end>m_numDim) throw myException("ContinuousProblem::initializeSolution@ContinuousProblem::initializePartSolution");
    CodeVReal&s=dynamic_cast< CodeVReal&>(result);
	double l,u;
	for(int i=begin;i<end;i++){
		m_searchRange.getSearchRange(l,u,i);
		s[i]= l+(u-l)*Global::msp_global->mp_uniformAlg->Next();
	}	
}
void ContinuousProblem::initializeSolution(VirtualEncoding &s,const int idx,const int maxIdx){
	CodeVReal&result=dynamic_cast< CodeVReal&>(s);
	switch(m_popInitialMode){
		case POP_INIT_UNIFORM:{
			//uniformly distributed in the entire search space
			double l,u;
			for(int i=0;i<m_numDim;i++){
				m_searchRange.getSearchRange(l,u,i);
				result[i]= l+(u-l)*Global::msp_global->mp_uniformAlg->Next();
			}
			}
			break;
		case POP_INIT_ORTHONORM:
			// orthonomalization method
			{///warning: need to be checked for verification
				double l,u;
				for( int i=0;i<m_numDim;i++){
					m_searchRange.getSearchRange(l,u,i);
					if(idx==0){
						result[i]= l;
					}else if(idx<maxIdx-1){
						result[i]= l+idx*((u-l)/(maxIdx-1));
					}else{
						result[i]= u;
					}
				}
			}
			break;
		case POP_INIT_CENTER:
			// nomarlly distributed in the center of the search space with variance of half domain
			{
				double l,u;
				for( int i=0;i<m_numDim;i++){
				m_searchRange.getSearchRange(l,u,i);
				result[i]= (l+u)/2+Global::msp_global->mp_normalAlg->NextNonStand(0,(u-l)/2);
			}

			}
			break;
		case POP_INIT_USER_DEFINED:
		///TO DO
			break;
	}
	
}
void ContinuousProblem::initializeSolution(const VirtualEncoding &ref,VirtualEncoding &s,double range){
	const CodeVReal&base=dynamic_cast< const CodeVReal&>(ref);
	CodeVReal &result=dynamic_cast< CodeVReal&>(s);
	double l,u;
	for( int i=0;i<m_numDim;i++){
		m_searchRange.getSearchRange(l,u,i);
		if(range==0) 		result[i]=(base[i])+Global::msp_global->mp_normalAlg->Next();
		else result[i]=(base[i])+Global::msp_global->mp_normalAlg->NextNonStand(0,range);

		if((result[i])<l)  result[i]=Global::msp_global->getRandFloat(l,(base[i]));
		if((result[i])>u) result[i]=Global::msp_global->getRandFloat((base[i]),u);
	}
}

ContinuousProblem::ContinuousProblem(const int rId, const int rDimNumber,  string rName, int numObj):\
	Problem(rId,rDimNumber,rName,numObj),m_disAccuracy(0.1),m_searchRange(rDimNumber),m_globalOpt(rDimNumber,numObj){		
		addProTag(CONT);
		if (m_id==Global::ms_curProId)		Solution<CodeVReal>::allocateMemoryWB(rDimNumber,numObj);		
}
ContinuousProblem& ContinuousProblem::operator=(const ContinuousProblem &rhs){
	if(this==&rhs) return *this;

	Problem::operator=(rhs);

	m_searchRange=rhs.m_searchRange;
	m_disAccuracy=rhs.m_disAccuracy;
	m_globalOpt=rhs.m_globalOpt;
	return *this;
}
void ContinuousProblem::parameterSetting(Problem * rP){
	Problem::parameterSetting(rP);
	
	//m_searchRange=dynamic_cast<ContinuousProblem*>(rP)->m_searchRange;
	m_disAccuracy=dynamic_cast<ContinuousProblem*>(rP)->m_disAccuracy;
	//m_globalOpt=dynamic_cast<ContinuousProblem*>(rP)->m_globalOpt;
}
bool ContinuousProblem::isSame(const VirtualEncoding &ss1, const VirtualEncoding & ss2){
	const CodeVReal &s1=dynamic_cast<const CodeVReal&>(ss1);
	const CodeVReal &s2=dynamic_cast<const CodeVReal&>(ss2);
	bool flag=true;
		for( int j=0;j<GET_NUM_DIM;j++){
		if((s1[j])!=(s2[j])&&fabs((double)(s1[j]-s2[j]))>1.E10-10){ 
			flag=false;
			break;
		}
		}
	return flag; 
	
}

ContinuousProblem::~ContinuousProblem(){

};

double ContinuousProblem::getDomainSize(){
	return m_searchRange.getDomainSize();
}
void ContinuousProblem::setDisAccuracy(double dis){
	m_disAccuracy=dis;
}


void ContinuousProblem::setSearchRange(double l, double u){
	m_searchRange.setBoundary(l,u);
}
void ContinuousProblem::setSearchRange(const vector<double> &l, const vector<double> &u){
	m_searchRange.setBoundary(l,u);
}
bool ContinuousProblem::getBoundaryFlag(int i){
	return m_searchRange.getBoundaryFlag(i);
}

double ContinuousProblem::getDistance(const VirtualEncoding &ss1, const VirtualEncoding &ss2, DistanceMode mode){
	
	const CodeVReal &s1=dynamic_cast<const CodeVReal&>(ss1);
	const CodeVReal &s2=dynamic_cast<const CodeVReal&>(ss2);
	double dis=0;
	
	switch(mode){
		case DIS_EUCLIDEAN:
		{
			for(int i=0;i<m_numDim;i++){
				dis+=(double)((s1[i]-s2[i])*(s1[i]-s2[i]));
			}
			return sqrt(dis);
		}
		case DIS_MANHATTAN:
		{
			for(int i=0;i<m_numDim;i++){
				dis+=fabs((double)(s1[i]-s2[i]));
			}
			return dis;
		}
		case DIS_HAMMING:
		{
			for(int i=0;i<m_numDim;i++){
				if(s1[i]!=s2[i]) dis+=1;	
			}
			return dis;
		}
		default: throw myException("undefined distance mode @Problem::getDistance");
	}

}

void ContinuousProblem::allocateMemory(const int numDim){
	m_globalOpt.resizeDim(numDim);
	m_searchRange.resizeDim(numDim);
}

bool ContinuousProblem::getObjGlobalOpt(vector<vector<double>> &value){
	if(m_globalOpt.flagGloObj()){
		value.clear();
		for(unsigned i=0;i<m_globalOpt.getNumOpt();i++)	value.push_back(m_globalOpt[i].obj());
        return true;
    }
    return false;

}
bool ContinuousProblem::getObjGlobalOpt(vector<double> &value){
	if(m_globalOpt.flagGloObj()){
		value=m_globalOpt[0].obj();
        return true;
    }
	return false;
	
}


bool ContinuousProblem::isGlobalOptKnown(){
	return m_globalOpt.flagLoc();
}

void ContinuousProblem::resizeDim(int num){
	Problem::resizeDim(num);
	m_searchRange.resizeDim(num);
	m_globalOpt.resizeDim(num);
}
void ContinuousProblem::resizeObj(int num){
	Problem::resizeObj(num);
	m_globalOpt.resizeObj(num);
}
