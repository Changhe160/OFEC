
// Generate Pareto front for DTLZ1-4 functions following NSGA-III's paper.
//
// An Evolutionary Many-Objective Optimization Algorithm Using Reference-point Based Non-dominated Sorting Approach, Part I:
// Solving Problems with Box Constraints
//
// http://dx.doi.org/10.1109/TEVC.2013.2281535

#include "DTLZ.h"

typedef vector<double> TObjVec;
typedef vector<TObjVec> TFront;

void GeneratePF_OneLayer(ostream &ofile, const string &problem_name, int M, int p);
void GeneratePF_TwoLayers(ostream &os, const string &problem_name, int M, int outside_p, int inside_p);

DTLZ::DTLZ(int ID, int numDim, const string &proName, int numObj) :BenchmarkFunction(ID, numDim, proName, numObj)
{
	if (m_numObj > m_numDim) throw myException("the number of dim must be greater or eaqual to the number of obj for DTLZ pros");
    setSearchRange(0.,1.);
	vector<ProTag> p_tag(1,MOP);
	p_tag.push_back(CONT);
	setProTag(p_tag);
	setOptType(MIN_OPT,-1);
	m_popInitialMode=POP_INIT_UNIFORM;
	generateAdLoadPF();
}

void DTLZ::generateAdLoadPF()
{
    const string problem_name[]= {"FUN_MOP_DTLZ1", "FUN_MOP_DTLZ2", "FUN_MOP_DTLZ3", "FUN_MOP_DTLZ4"};
	stringstream os;
	os<<Global::g_arg[param_workingDir]<<"Problem/FunctionOpt/Data/PF_"<<Global::g_arg[param_proName]<<"("<<Global::g_arg[param_numObj]<<")"<<"_Opt.txt";

    for (int i=0; i<4; i+=1) // problem
    {
		if(Global::g_arg[param_proName]!=problem_name[i])
			continue;
        const int M[5] = {3, 5, 8, 10, 15};
        for (int j=0; j<5; j+=1) // objectives
        {
			if(Global::g_arg[param_numObj]!=M[j])
				continue;
			ifstream infile(os.str());
			if(infile) 
			{	infile.close(); break; }
            ofstream ofile(os.str());

            if (M[j] <= 5) // #objectives <= 5
            {
                int p[2] = {12, 6}; // Check Section V, Table I in the original paper
                GeneratePF_OneLayer(ofile, problem_name[i], M[j], p[j]);
            }
            else
            {
                int p[3][2] = {{3, 2}, {3, 2}, {2, 1}}; // Check Section V, Table I in the original paper
                GeneratePF_TwoLayers(ofile, problem_name[i], M[j], p[j-2][0], p[j-2][1]);
            }
			ofile.close();
        }
    }
	int numObj=Global::g_arg[param_numObj];
	ifstream infile(os.str());
	if(!infile)
		throw myException("please set your own pareto front @DTLZ::generatePF()");
	string str;
	int line=0;
	while(getline(infile,str))
		++line;
	m_globalOpt.setNumOpts(line);
	m_originalGlobalOpt.setNumOpts(line);
	m_originalGlobalOpt.setFlagLocTrue();
	infile.close();
	infile.clear();
	infile.open(os.str());
	for(int i=0;i<line;i++)
		for(int j=0;j<numObj;j++)
			infile>>m_originalGlobalOpt[i].data().m_obj[j];
	m_globalOpt=m_originalGlobalOpt;
	infile.close();
}

// ----------------------------------------------------------------------
void generate_recursive(TFront *pf, TObjVec *pt, size_t num_objs,
                        size_t left, size_t total, size_t element)
{
    if (element == num_objs-1)
    {
        (*pt)[element] = left;
        pf->push_back(*pt);
    }
    else
    {
        for (size_t i=0; i<=left; i+=1)
        {
            (*pt)[element] = i;
            generate_recursive(pf, pt, num_objs, left-i, total, element+1);
        }
    }
}
// ----------------------------------------------------------------------
void GenerateWeight(TFront *pf, size_t M, size_t p)
{
    TObjVec pt(M);

    generate_recursive(pf, &pt, M, p, p, 0);
}
// ----------------------------------------------------------------------
void GeneratePF_OneLayer(ostream &os, const string &problem_name, int M, int p)
{
    TFront PF;

    int num_objectives = M, num_divisions = p;
    GenerateWeight(&PF, num_objectives, num_divisions);

    if (problem_name == "FUN_MOP_DTLZ1")
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << (0.5*PF[i][j])/num_divisions << ' ';
            }
            os << endl;
        }
    }
    else // DTLZ2-4
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            double sum = 0;

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                sum += PF[i][j]*PF[i][j];
            }

            double k = sqrt(1.0/sum);

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << k*PF[i][j] << ' ';
            }
            os << endl;
        }
    } // else

}// GeneratePF_OneLayer()
// ----------------------------------------------------------------------
void GeneratePF_TwoLayers(ostream &os, const string &problem_name, int M, int outside_p, int inside_p)
{

    GeneratePF_OneLayer(os, problem_name, M, outside_p);


    TFront PF;

    int num_objectives = M, num_divisions = inside_p;
    GenerateWeight(&PF, num_objectives, num_divisions);

    for (size_t i=0; i<PF.size(); i+=1)
    {
        for (size_t j=0; j<PF[i].size(); j+=1)
        {
            PF[i][j] = (static_cast<double>(num_divisions)/M+PF[i][j])/2; // (k=num_divisions/M, k, k, ..., k) is the center point
        }
    }

    if (problem_name == "FUN_MOP_DTLZ1")
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << (0.5*PF[i][j])/num_divisions << ' ';
            }
            os << endl;
        }
    }
    else // DTLZ2-4
    {
        for (size_t i=0; i<PF.size(); i+=1)
        {
            double sum = 0;

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                sum += PF[i][j]*PF[i][j];
            }

            double k = sqrt(1.0/sum);

            for (size_t j=0; j<PF[i].size(); j+=1)
            {
                os << k*PF[i][j] << ' ';
            }
            os << endl;
        }
    } // else

}// GeneratePF_TwoLayers()
