#include "LKH.h"
#include "Genetic.h"

/*
 * This file contains the main function of the program.
 */
vector<double> LKH::LKHAlg::mv_cost;
LKH::LKHAlg::LKHAlg(ParamMap &v) :Algorithm(-1,string())
{
	initializeParam();
	string file=v[param_dataFile1];
	ProblemFileName=new char[file.size()+1];
	ParameterFileName=new char[file.size()+1];
	int loc=file.find(".tsp");
	file.erase(loc,4);
	strcpy(ProblemFileName,file.c_str());
	strcpy(ParameterFileName,ProblemFileName);
	mv_cost.resize(int(MAX_NUM_RUN));
	/* Read the specification of the problem */
    ReadParameters();
    MaxMatrixDimension = 10000;
    ReadProblem();
}

LKH::LKHAlg::LKHAlg(const string &fileName) :Algorithm(-1,string())
{
	initializeParam();
	ProblemFileName=new char[fileName.size()+1];
	ParameterFileName=new char[fileName.size()+1];
	strcpy(ProblemFileName,fileName.c_str());
	strcpy(ParameterFileName,fileName.c_str());
	mv_cost.resize(int(MAX_NUM_RUN));
	/* Read the specification of the problem */
    ReadParameters();
    MaxMatrixDimension = 10000;
    ReadProblem();
}

LKH::LKHAlg::~LKHAlg()
{
	delete []ProblemFileName;
	delete []ParameterFileName;
}

void LKH::LKHAlg::initializeParam()
{
	savePtr=0;
	BestTour=0;  /* Table containing best tour found */
	BetterTour=0;        /* Table containing the currently best tour in a run */
	CacheVal=0;  /* Table of cached distances */
	CacheSig=0;  /* Table of the signatures of cached distances */
	CostMatrix=0;        /* Cost matrix */
	FirstActive=0; LastActive=0; /* First and last node in the list of "active" nodes */
	FirstNode=0;        /* First node in the list of nodes */
	FirstSegment=0;  /* A pointer to the first segment in the cyclic list of segments */
	FirstSSegment=0;        /* A pointer to the first super segment in the cyclic list of segments */
	Heap=0;    /* Heap used for computing minimum spanning trees */
	HTable=0;      /* Hash table used for storing tours */
	LastLine=0; /* Last input line */
	NodeSet=0;  /* Array of all nodes */
	Rand=0; /* Table of random values */
	SwapStack=0;  /* Stack of SwapRecords */

/* The following variables are read by the functions ReadParameters and 
   ReadProblem: */

    ParameterFileName=0; ProblemFileName=0; PiFileName=0; TourFileName=0; OutputTourFileName=0; InputTourFileName=0;
    CandidateFileName=0; InitialTourFileName=0; SubproblemTourFileName=0; MergeTourFileName=0; Name=0; Type=0; EdgeWeightType=0; 
	EdgeWeightFormat=0; EdgeDataFormat=0; NodeCoordType=0; DisplayDataType=0;

    ParameterFile=0; ProblemFile=0; PiFile=0; InputTourFile=0; TourFile=0; InitialTourFile=0; SubproblemTourFile=0; MergeTourFile=0;
    Distance=0; D=0; C=0; c=0; BestMove=0; BacktrackMove=0; BestSubsequentMove=0;
}

ReturnFlag LKH::LKHAlg::run_()
{
	GainType Cost, OldOptimum;
    double Time, LastTime = GetTime();

    if (SubproblemSize > 0) {
        if (DelaunayPartitioning)
            SolveDelaunaySubproblems();
        else if (KarpPartitioning)
            SolveKarpSubproblems();
        else if (KCenterPartitioning)
            SolveKCenterSubproblems();
        else if (KMeansPartitioning)
            SolveKMeansSubproblems();
        else if (RohePartitioning)
            SolveRoheSubproblems();
        else if (MoorePartitioning || SierpinskiPartitioning)
            SolveSFCSubproblems();
        else
            SolveTourSegmentSubproblems();
    }
    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType) LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        WriteTour(OutputTourFileName, BestTour, BestCost);
        WriteTour(TourFileName, BestTour, BestCost);
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
        if (*MaxPopulationSize > 1) {
            /* Genetic algorithm */
            int i;
            for (i = 0; i < *PopulationSize; i++) {
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i,this);
                if (TraceLevel >= 1 && Cost < OldCost) {
                    printff("  Merged with %d: Cost = " GainFormat, i + 1,
                            Cost);
                    if (Optimum != MINUS_INFINITY && Optimum != 0)
                        printff(", Gap = %0.4f%%",
                                100.0 * (Cost - Optimum) / Optimum);
                    printff("\n");
                }
            }
            if (!HasFitness(Cost)) {
                if (*PopulationSize < *MaxPopulationSize) {
                    AddToPopulation(Cost,this);
                    if (TraceLevel >= 1)
                        PrintPopulation(this);
                } else if (Cost < Fitness.get()[*PopulationSize - 1]) {
                    i = ReplacementIndividual(Cost,this);
                    ReplaceIndividualWithTour(i, Cost,this);
                    if (TraceLevel >= 1)
                        PrintPopulation(this);
                }
            }
        } else if (Run > 1)
            Cost = MergeBetterTourWithBestTour();
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
            WriteTour(OutputTourFileName, BestTour, BestCost);
            WriteTour(TourFileName, BestTour, BestCost);
        }
        OldOptimum = Optimum;
        if (Cost < Optimum) {
            if (FirstNode->InputSuc) {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
			printff("Run %d: Cost = " GainFormat, Global::msp_global->m_runId, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum);
        //    printff(", Time = %0.2f sec. %s\n\n", Time,
         //           Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
			printff(", Time = %0.2f sec. \n", Time);
        }
        if (StopAtOptimum && Cost == OldOptimum && *MaxPopulationSize >= 1) {
            Runs = Run;
            break;
        }
        if (*PopulationSize >= 2 &&
            (*PopulationSize == *MaxPopulationSize ||
             Run >= 2 * *MaxPopulationSize) && Run < Runs) {
            Node *N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(*PopulationSize, 1.25,this);
            do
                Parent2 = LinearSelection(*PopulationSize, 1.25,this);
            while (Parent2 == Parent1);
            ApplyCrossover(Parent1, Parent2,this);
            N = FirstNode;
            do {
                int d = (this->*C)(N, N->Suc);
                AddCandidate(N, N->Suc, d, INT_MAX);
                AddCandidate(N->Suc, N, d, INT_MAX);
                N = N->InitialSuc = N->Suc;
            }
            while (N != FirstNode);
        }
		mv_cost[Global::msp_global->m_runId]=Cost;
        SRandom(++Seed);
    }
 //   PrintStatistics();
	freeAll();
	FreeStructures();
	return Return_Terminate;
}

void LKH::LKHAlg::freeAll()
{
	freeERXT();
	freeGain23();
	freeSequence();
	freeCreateDelaunay();
	freeCreateQuadrant();
	freeGreedyTour();
	freePatchCycle();
	freeReadPenalties();
}
