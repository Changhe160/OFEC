/*************************************************************************
* Project: Library of Evolutionary Algoriths
*************************************************************************
* Author: Changhe Li & Ming Yang & Yong Xia
* Email: changhe.lw@google.com Or yangming0702@gmail.com Or cugxiayong@gmail.com
* Language: C++
*************************************************************************
*  This file is part of EAlib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Oct 2014
// Last modified:


#include "TravellingSalesman.h"
#include "../../../Global/global.h"
#include "../../../Algorithm/Other/LKH/LKH.h"
#include<string>

using namespace std;

static boost::mutex g_mutex1;

TravellingSalesman::TravellingSalesman(ParamMap &v):Problem((v[param_proId]), (v[param_numDim]),(v[param_proName]),1),mvv_solut(m_numDim,vector<bool>(m_numDim,false)),m_globalOpt(m_numDim,m_numObj) 
{
	v[param_numObj]=1;
	mvvv_cost.resize(m_numObj);
	for(int i=0;i<m_numObj;i++)
		mvvv_cost[i].resize(m_numDim);
	for(int i=0;i<m_numObj;i++)
		for(int j=0;j<m_numDim;j++)
			mvvv_cost[i][j].resize(m_numDim);
	m_globalOpt.setFlagLocTrue();
	string file=v[param_proFileName];
	m_fileName=file;
	int loc=m_fileName.find(".tsp");
	m_fileName.erase(loc,4);
	readProblem();

	m_searchRange.setBoundary(0,(v[param_numDim])-1);
	int mode=int(v[param_populationInitialMethod]);
	m_popInitialMode=static_cast<PopInitMethod>(mode);
	addProTag(TSP);
	setOptType(MIN_OPT);
	
	v[param_sampleFre]=m_numDim*2;
	
	//g_mutex1.lock();
	Solution<CodeVInt>::allocateMemoryWB(m_numDim,m_numObj);
	//g_mutex1.unlock();
	
}

TravellingSalesman::~TravellingSalesman()
{
}

TravellingSalesman::TravellingSalesman(const int rId, const int rDimNumber, string rName, string fileName, int numObj):Problem(rId,rDimNumber,rName,numObj),mvv_solut(m_numDim,vector<bool>(m_numDim,false)),m_globalOpt(m_numDim,m_numObj)  
{
	
	mvvv_cost.resize(m_numObj);
	for(int i=0;i<m_numObj;i++)
		mvvv_cost[i].resize(m_numDim);
	for(int i=0;i<m_numObj;i++)
		for(int j=0;j<m_numDim;j++)
			mvvv_cost[i][j].resize(m_numDim);
	m_globalOpt.setFlagLocTrue();
	readProblem();

	m_searchRange.setBoundary(0,rDimNumber-1);
	int mode=Global::g_arg[param_populationInitialMethod];
	m_popInitialMode=static_cast<PopInitMethod>(mode);
	addProTag(TSP);
	setOptType(MIN_OPT);
	g_mutex1.lock();
	Solution<CodeVInt>::allocateMemoryWB(m_numDim,m_numObj);
	g_mutex1.unlock();
}


void TravellingSalesman::readProblem()
{
	size_t i;
	string Line;
	char *edgeType=0,*edgeFormat=0;
	char *Keyword=0;
	char *Type=0,*Format=0;
	const char *Delimiters=" ():=\n\t\r\f\v\xef\xbb\xbf";
	vector<vector<double> > coordinate;
	ostringstream oss1,oss2;
	ifstream infile;
	oss1<<Global::g_arg[param_workingDir]<<"Problem/Combination/TSP/data/"<<m_fileName<<".tsp";
	infile.open(oss1.str().c_str());
	if(!infile){
		throw myException("read travelling salesman data error");
	}
	char *savePtr;
	while(getline(infile,Line))
	{
		if(!(Keyword=gStrtok_r((char*)Line.c_str(),Delimiters,&savePtr)))
			continue;
		for(i=0;i<strlen(Keyword);i++)
			Keyword[i]=toupper(Keyword[i]);
		if(!strcmp(Keyword,"NAME"))
			continue;
		else if(!strcmp(Keyword,"COMMENT"))
			continue;
		else if(!strcmp(Keyword,"TYPE"))
			continue;
		else if(!strcmp(Keyword,"DIMENSION"))
		{
			char *token = gStrtok_r(0,Delimiters,&savePtr);
			m_numDim=atoi(token);
		}
		else if(!strcmp(Keyword,"EDGE_WEIGHT_TYPE"))
		{
			edgeType=gStrtok_r(0,Delimiters,&savePtr);
			Type=new char[strlen(edgeType)+1];
			strcpy(Type,edgeType);
		}
		else if(!strcmp(Keyword,"EDGE_WEIGHT_FORMAT"))
		{
			edgeFormat=gStrtok_r(0,Delimiters,&savePtr);
			Format=new char[strlen(edgeFormat)+1];
			strcpy(Format,edgeFormat);
		}
		else if(!strcmp(Keyword,"NODE_COORD_SECTION"))
		{
			i=0;
			vector<double> temp(2);
			while(infile>>Line)
			{
				infile>>temp[0];
				infile>>temp[1];
				coordinate.push_back(temp);
				++i;
				if(i==(size_t)m_numDim) break;
			}
			break;
		}
		else if(!strcmp(Keyword,"EDGE_WEIGHT_SECTION"))
		{
			if(!strcmp(Format,"LOWER_DIAG_ROW"))
			{
				for(int i=0;i<m_numDim;i++)
					for(int j=0;j<=i;j++)
					{
						infile>>mvvv_cost[0][i][j];
						mvvv_cost[0][j][i]=mvvv_cost[0][i][j];
					}
			}
			else if(!strcmp(Format,"UPPER_DIAG_ROW"))
			{
				for(int i=0;i<m_numDim;i++)
					for(int j=i;j<m_numDim;j++)
					{
						infile>>mvvv_cost[0][i][j];
						mvvv_cost[0][j][i]=mvvv_cost[0][i][j];
					}
			}
			else if(!strcmp(Format,"FULL_MATRIX"))  
			{
				for(int i=0;i<m_numDim;i++)
					for(int j=0;j<m_numDim;j++)
						infile>>mvvv_cost[0][i][j];
			}
			else 
				throw myException("no exiting this edgeFormat in function readProblem in TravellingSalesman.cpp, please add it here!");
		}
	}
	infile.close();
	infile.clear();
	if(!Type)
		throw myException("file format error in function readProblem in TravellingSalesman.cpp");
	if(strcmp(Type,"EXPLICIT"))
		calculateEdgeWeight(Type,coordinate);
	if(Type) delete []Type;
	if(Format) delete []Format;
	oss2<<Global::g_arg[param_workingDir]<<"Problem/Combination/TSP/data/"<<m_fileName<<".opt.tour";
	infile.open(oss2.str().c_str());
	if(!infile)
		throw myException("read travelling salesman optimal tour file error");
	string oldLine;
	while(getline(infile,Line))
	{
		oldLine=Line;
		if(!(Keyword=gStrtok_r((char*)Line.c_str(),Delimiters,&savePtr)))
			continue;
		for(i=0;i<strlen(Keyword);i++)
			Keyword[i]=toupper(Keyword[i]);
		if(!strcmp(Keyword,"TOUR_SECTION"))
		{
			int temp;
			for(i=0;i<(size_t)m_numDim;i++)
			{
				infile>>temp;
				--temp;
				m_globalOpt[0].data().m_x[i]=TypeVar(temp);
			}
			break;
		}
		else if(!strcmp(Keyword,"COMMENT"))  //get optimal cost
		{
			size_t j=0;
			for(int z=oldLine.size()-1;z>=0;z--)
				Line[j++]=oldLine[z];
			Keyword=gStrtok_r((char*)Line.c_str(),Delimiters,&savePtr);
			size_t len=strlen(Keyword)-1;
			for(i=0;i<=len;i++)
				oldLine[i]=Keyword[len-i];
			oldLine[len+1]='\0';
			vector<vector<double> > obj;
			obj.push_back(vector<double>(1,atof(oldLine.c_str())));
			m_globalOpt.setGloObj(obj);
			m_globalOpt.setFlagLocTrue();
		}
		else 
			continue;
	}
	infile.close();
	infile.clear();
}

void TravellingSalesman::calculateEdgeWeight(char *edgeType,vector<vector<double> >& coordinate)
{
	if(!strcmp(edgeType,"EUC_2D"))
	{
		for(int i=0;i<m_numDim;i++)
			for(int j=0;j<m_numDim;j++)
				mvvv_cost[0][i][j]=static_cast<int>(sqrt(pow(coordinate[i][0]-coordinate[j][0],2)+pow(coordinate[i][1]-coordinate[j][1],2))+0.5); 
	}
	else if(!strcmp(edgeType,"ATT"))
	{
		double r;
		int t;
		for(int i=0;i<m_numDim;i++)
		{
			for(int j=0;j<m_numDim;j++)
			{
				r=sqrt((pow(coordinate[i][0]-coordinate[j][0],2)+pow(coordinate[i][1]-coordinate[j][1],2))/10.0);
				t=static_cast<int>(r+0.5);
				if(t<r) mvvv_cost[0][i][j]=t+1;
				else mvvv_cost[0][i][j]=t;
			}
		}
	}
	else if(!strcmp(edgeType,"GEO"))
	{
		double pi=3.141592,RRR=6378.388,min;
		int deg;
		for(int i=0;i<m_numDim;i++)
		{
			for(int j=0;j<2;j++)
			{
				deg=static_cast<int>(coordinate[i][j]);
				min=coordinate[i][j]-deg;
				coordinate[i][j]=pi*(deg+5.0*min/3.0)/180.0;
			}
		}
		double q1,q2,q3;
		for(int i=0;i<m_numDim;i++)
		{
			for(int j=0;j<m_numDim;j++)
			{
				q1=cos(coordinate[i][1]-coordinate[j][1]);
				q2=cos(coordinate[i][0]-coordinate[j][0]);
				q3=cos(coordinate[i][0]+coordinate[j][0]);
				mvvv_cost[0][i][j]=static_cast<int>(RRR*acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0);
			}
		}
	}
	else
		throw myException("no exiting this edgeType in function calculateEdgeWeight in TravellingSalesman.cpp, please add it here!");
}

ReturnFlag TravellingSalesman::evaluate_(VirtualEncoding &s_, bool rFlag, ProgramMode mode,  bool flag2)
{	
	CodeVInt &s=dynamic_cast< CodeVInt&>(s_);

	for(int i=0;i<m_numObj;i++)
		s.m_obj[i]=0;
	int row,col;
	for(int n=0;n<m_numObj;n++)
	{
		for(size_t i=0;i<m_numDim;i++)
		{
			row=s.m_x[i];
			col=s.m_x[(i+1)%m_numDim];
			s.m_obj[n]+=mvvv_cost[n][row][col];
		}
	}
	if (flag2){
		if (rFlag)	m_evals++;

		if (Global::msp_global->mp_algorithm.get() != nullptr&&Global::msp_global->mp_algorithm->ifTerminating()) return Return_Terminate;
		return Return_Normal;
	}
	return Return_Normal;
}

bool TravellingSalesman::isValid(const VirtualEncoding &s1)
{
	const CodeVInt& s=dynamic_cast<const CodeVInt&>(s1);

	for(int i=0;i<m_numDim;i++){  //judge the range		
		if((s[i])<m_searchRange.m_lower||(s[i])>m_searchRange.m_upper) return false;		
	}
	vector<int> flag(m_numDim,0);  //judge whether has the same gene
	int temp;
	for(int i=0;i<m_numDim;i++)
	{
		temp=s[i];
		flag[temp]=1;
	}
	for(int i=0;i<m_numDim;i++)
		if(flag[i]==0)
			return false;
	return true;
}


bool TravellingSalesman::isSame(const VirtualEncoding &ss1, const VirtualEncoding &ss2)
{
	const CodeVInt& s1=dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2=dynamic_cast<const CodeVInt&>(ss2);
	int i,j,pos;
	for(i=0;i<m_numDim;i++)
	{
		if(s1[0]==s2[i])
		{
			pos=i;
			break;
		}
	}
	j=pos;
	for(i=0;i<m_numDim;i++)
	{
		if(s1[i]!=s2[j])  break;
		else j=(j+1)%m_numDim;
	}
	if(i==m_numDim) return true;
	j=pos;
	for(i=0;i<m_numDim;i++)
	{
		if(s1[i]!=s2[j])  break;
		else j=(j-1+m_numDim)%m_numDim;
	}
	if(i==m_numDim) return true;
	return false;
}

void TravellingSalesman::initializeSolution(VirtualEncoding &result_,const int idx,const int maxId)
{
	CodeVInt& result=dynamic_cast< CodeVInt&>(result_);

	static vector< vector<int> > candidateSet(0);       //candidate set
	static vector< vector<int> >  nearby(0);         //nearby cities of some city
	static string file("");
	g_mutex1.lock();
	double pg,pr;
	if(Global::g_arg.find(param_interTest1)!=Global::g_arg.end()) pg=Global::g_arg[param_interTest1];
	else pg=0.5;

	if(Global::g_arg.find(param_interTest2)!=Global::g_arg.end()) pr=Global::g_arg[param_interTest2];
	else pr=0.2;

	if(file!=m_fileName)
	{
		file=m_fileName;
		candidateSet.resize(0);
		nearby.resize(0);
	}
	g_mutex1.unlock();
	switch(m_popInitialMode){
	case POP_INIT_HEURIS:
	{
		//utilize candidate set, top and randomness to initialize a solution
		g_mutex1.lock();
		if(!candidateSet.size())
		{
			nearby.resize(m_numDim);
			candidateSet.resize(m_numDim);
			for(int i=0;i<m_numDim;i++)
				nearby[i].resize(int(0.15*m_numDim));
			findNearbyCity(nearby);
			createCandidateSets(candidateSet);
		}
		g_mutex1.unlock();
		int i,j;
		int pos;
		int num=nearby[0].size();
		vector<int> visited(m_numDim,0);
		result[0]=int((m_numDim-1)*Global::msp_global->mp_uniformAlg->Next());  
		visited[int(result[0])]=1;
		for(i=1;i<m_numDim;i++)
		{
			double p=Global::msp_global->mp_uniformAlg->Next();
			if(p<pg) //TOP, Greedy strategy  0.5
			{
				j=0;
				pos=result[i-1];
				while(j<num)
				{
					if(visited[nearby[pos][j]]==1)
							j++;
					else 
					{
						result[i]=nearby[pos][j];
						visited[nearby[pos][j]]=1;
						break;
					}
				} 
				if(j==num) //如果遍历所有num个城市都不行，随机给一个城市
					result[i]=selectCityRandom(visited,m_numDim);
			}
			else if(p<pg+pr)  //TOP, select one city from 15% nearest cities randomly 0.7
			{
				j=selectCityRandom(nearby,visited,num,result[i-1]);
				if(j!=-1) result[i]=j;
				else //如果遍历所有num个城市都不行，随机给一个城市
					result[i]=selectCityRandom(visited,m_numDim);
			}
			else  //candidate set Strategy
			{
				j=0;
				pos=result[i-1];
				while(j<candidateSet[pos].size())
				{
					if(visited[candidateSet[pos][j]]==1)
							j++;
					else
					{
						result[i]=candidateSet[pos][j];
						visited[candidateSet[pos][j]]=1;
						break;
					}
				}
				if(j==candidateSet[pos].size())  //候选集结束还没找到下个点
					result[i]= selectCityRandom(visited,m_numDim);
			}
		}
	}
	break;

	case POP_INIT_RANDOM:
		{
			vector<int> temp;
			int i,pos,num=result.m_x.size();
			for(i=0;i<num;i++)
				temp.push_back(int(i));
			for(i=0;i<num;i++)
			{
				pos=int((num-1-i)*Global::msp_global->mp_uniformAlg->Next());
				result[i]=temp[pos];
				temp[pos]=temp[num-1-i];
			}
		}
		break;
	}
	if(!isValid(result))
		throw myException("error in @TravellingSalesman::initializeSolution() in TravellingSalesman.cpp");
}

int selectCityRandom(vector< vector<int> >& matrix,vector<int> &visited,int num,int row)
{
	int i,pos,flag=-1;
	vector<int> arr(num);
	for(i=0;i<num;i++)
		arr[i]=matrix[row][i];
	i=num-1;
	while(i>=0)
	{
		pos=int(i*Global::msp_global->mp_uniformAlg->Next());
		if(visited[arr[pos]]==1)
		{
			arr[pos]=arr[i];
			i--;
		}
		else 
		{
			visited[arr[pos]]=1;
			flag=arr[pos];
			break;
		}
	}
	return flag;
}

int selectCityRandom(vector<int> &visited,int dim)
{
	int i,pos,flag=-1;
	vector<int> arr(dim);
	for(i=0;i<dim;i++)
		arr[i]=i;
	i=dim-1;
	while(i>=0)
	{
		pos=int(i*Global::msp_global->mp_uniformAlg->Next());
		if(visited[arr[pos]]==1)
		{
			arr[pos]=arr[i];
			i--;
		}
		else 
		{
			visited[arr[pos]]=1;
			flag=arr[pos];
			break;
		}
	}
	return flag;
}

void TravellingSalesman::findNearbyCity(vector<vector<int> > &nearby,int n)
{
	int i,j,z,pos;
	double min;
	vector<int> visited(m_numDim);
	int num=nearby[0].size();
	for(i=0;i<m_numDim;i++)
	{
		for(j=0;j<m_numDim;j++)
			visited[j]=0;
		for(j=0;j<num;j++)
		{
			min=0xfffffff;
			for(z=0;z<m_numDim;z++)
			{
				if(min>mvvv_cost[n][i][z]&&i!=z&&visited[z]==0)
				{
					min=mvvv_cost[n][i][z];
					pos=z;
				}
			}
			visited[pos]=1;
			nearby[i][j]=pos;
		}
	}
}

void TravellingSalesman::createCandidateSets(vector<vector<int> > &candidateSets)
{
	LKH::LKHAlg lkh(m_fileName);
	lkh.AllocateStructures();
    lkh.CreateCandidateSet();
	LKH::LKHAlg::Node * start=lkh.FirstNode;
	do
	{
		for(LKH::LKHAlg::Candidate *NFrom=start->CandidateSet; NFrom->To; NFrom++)
			candidateSets[start->Id-1].push_back(NFrom->To->Id-1);
	}while((start=start->Suc)!=lkh.FirstNode);
}


void TravellingSalesman::prim(vector<vector<int> > &mstEdge,int n)
{
	vector<int> near(m_numDim);
	int i,j,k,l;
	int Min=mvvv_cost[n][0][1];
	k=0;l=1;
	for(i=0;i<m_numDim;i++)  //find the nearest edge in the graph
	{
		for(j=0;j<m_numDim;j++)
		{
			if(Min>mvvv_cost[n][i][j]&&i!=j)
			{
				Min=mvvv_cost[n][i][j];
				k=i;l=j;
			}
		}
	}
	mstEdge[0][0]=k; mstEdge[0][1]=l;
	for(i=0;i<m_numDim;i++)  //initial near
	{
		if(mvvv_cost[n][i][l]<mvvv_cost[n][i][k])
			near[i]=l;
		else near[i]=k;
	}
	near[l]=-1; near[k]=-1;
	for(i=2;i<m_numDim;i++)//find the remain n-2 edges
	{
		Min=0xfffffff;
		for(int z=0;z<m_numDim;z++)
		{
			if(near[z]!=-1&&Min>mvvv_cost[n][z][near[z]])
			{
				Min=mvvv_cost[n][z][near[z]];
				j=z;
			}
		}
		mstEdge[i-1][0]=j; mstEdge[i-1][1]=near[j];
		near[j]=-1;
		for(k=0;k<m_numDim;k++)
			if(near[k]!=-1&&mvvv_cost[n][k][near[k]]>mvvv_cost[n][k][j])
				near[k]=j;
	}
}
double TravellingSalesman::getDistance(const VirtualEncoding &ss1, const VirtualEncoding &ss2, DistanceMode mode){
	const CodeVInt& s1=dynamic_cast<const CodeVInt&>(ss1);
	const CodeVInt& s2=dynamic_cast<const CodeVInt&>(ss2);
	for(int i=0;i<m_numDim;i++){
		for(int j=0;j<m_numDim;j++){
			mvv_solut[i][j]=false;
		}
	}
	for(int i=0;i<m_numDim;i++){
		mvv_solut[s1[i]][s1[(i+1)%m_numDim]]=true;
	}
	double dis=0;
	for(int i=0;i<m_numDim;i++){
		if(!mvv_solut[s2[i]][s2[(i+1)%m_numDim]])dis+=1;
	}
	return dis;
}

pair<int,int> TravellingSalesman::getNextCity(const VirtualEncoding &s_,int n){
	const CodeVInt& s=dynamic_cast<const CodeVInt&>(s_);
	for(int i=0;i<m_numDim;i++){
		if(s[i]==n) return pair<int,int>(s[(i-1+m_numDim)%m_numDim],s[(i+1)%m_numDim]);
	}
	throw myException("error@TravellingSalesman::getNextCity(VirtualEncoding &s,int n)");
} 

bool TravellingSalesman::getObjGlobalOpt(vector<vector<double>> &value){
	if(m_globalOpt.flagGloObj()){
		value.clear();
		for(unsigned i=0;i<m_globalOpt.getNumOpt();i++)	value.push_back(m_globalOpt[i].obj());
        return true;
    }
    return false;

}
bool TravellingSalesman::getObjGlobalOpt(vector<double> &value){
	if(m_globalOpt.flagGloObj()){
		value=m_globalOpt[0].obj();
        return true;
    }
	return false;
	
}

const Optima<CodeVInt> & TravellingSalesman::getGOpt()const{
	return m_globalOpt;
}

Optima<CodeVInt> & TravellingSalesman::getGOpt(){
	return m_globalOpt;
}

TravellingSalesman *TravellingSalesman::getTypePtr(){
	return this;
}
TravellingSalesman &TravellingSalesman::getTypeRef(){
	return *this;
}
bool TravellingSalesman::isGlobalOptKnown(){
	return m_globalOpt.flagLoc();
}