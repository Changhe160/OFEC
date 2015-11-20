#include "LKH.h"
#include "GeoConversion.h"

/*
 * The SolveKarpSubproblems function attempts to improve a given tour 
 * by means of Karp's partitioning scheme. 
 *
 * The overall region containing the nodes is subdivided into rectangles
 * with SubproblemSize nodes in each rectangle. Each rectangle together 
 * with the given tour induces a subproblem consisting of all nodes inside 
 * the rectangle, and with edges fixed between nodes that are connected 
 * by tour segments whose interior points are outside the rectangle.  
 *  
 * If an improvement is found, the new tour is written to TourFile. 
 * The original tour is given by the SubproblemSuc references of the nodes.
 */

static void KarpPartition(int start, int end, LKH::LKHAlg * Alg);
static void CalculateSubproblems(int start, int end, LKH::LKHAlg * Alg);

static boost::thread_specific_ptr<LKH::LKHAlg::Node *> KDTree;
static boost::thread_specific_ptr<GainType> GlobalBestCost, OldGlobalBestCost;
static boost::thread_specific_ptr<int> CurrentSubproblem, Subproblems;

void LKH::LKHAlg::SolveKarpSubproblems()
{
    Node *N;
    double EntryTime = GetTime();
	if(!GlobalBestCost.get())
	{
		GlobalBestCost.reset(new GainType(0));
		OldGlobalBestCost.reset(new GainType(0));
		CurrentSubproblem.reset(new int(0));
		Subproblems.reset(new int(0));
	}

    AllocateStructures();
    ReadPenalties();

    /* Compute upper bound for the original problem */
    *GlobalBestCost = 0;
    N = FirstNode;
    do {
        if (!Fixed(N, N->SubproblemSuc))
            *GlobalBestCost += (this->*Distance)(N, N->SubproblemSuc);
        N->Subproblem = 0;
    }
    while ((N = N->SubproblemSuc) != FirstNode);
    if (TraceLevel >= 1) {
        if (TraceLevel >= 2)
            printff("\n");
        printff("*** Karp partitioning *** [Cost = " GainFormat "]\n",
                *GlobalBestCost);
    }
    if (WeightType == GEO || WeightType == GEOM || 
        WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS) {
        N = FirstNode;
        do {
            N->Xc = N->X;
            N->Yc = N->Y;
            N->Zc = N->Z;
            if (WeightType == GEO || WeightType == GEO_MEEUS)
                GEO2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
            else
                GEOM2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
        } while ((N = N->SubproblemSuc) != FirstNode);
        CoordType = THREED_COORDS;
    }
	KDTree.reset(BuildKDTree(SubproblemSize));
    if (WeightType == GEO || WeightType == GEOM ||
        WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS) {
        N = FirstNode;
        do {
            N->X = N->Xc;
            N->Y = N->Yc;
            N->Z = N->Zc;
        } while ((N = N->SubproblemSuc) != FirstNode);
        CoordType = TWOD_COORDS;
    }

    *Subproblems = 0;
    CalculateSubproblems(0, Dimension - 1,this);
    *CurrentSubproblem = 0;
    KarpPartition(0, Dimension - 1,this);
    free(KDTree.get());
    printff("\nCost = " GainFormat, *GlobalBestCost);
    if (Optimum != MINUS_INFINITY && Optimum != 0)
        printff(", Gap = %0.4f%%",
                100.0 * (*GlobalBestCost - Optimum) / Optimum);
    printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
            *GlobalBestCost < Optimum ? "<" : *GlobalBestCost ==
            Optimum ? "=" : "");
    if (SubproblemBorders && *Subproblems > 1)
        SolveSubproblemBorderProblems(*Subproblems, &*GlobalBestCost);
}

/*
 * The KarpPartition function subidivides the overall region into 
 * rectangles and attempts to solve the induced subproblems. 
 */

static void KarpPartition(int start, int end, LKH::LKHAlg * Alg)
{
    if (end - start + 1 <= Alg->SubproblemSize) {
        int i;
        (*CurrentSubproblem)++;
        for (i = start; i <= end; i++)
            KDTree.get()[i]->Subproblem = *CurrentSubproblem;
        *OldGlobalBestCost = *GlobalBestCost;
        Alg->SolveSubproblem(*CurrentSubproblem, *Subproblems, &*GlobalBestCost);
        if (Alg->SubproblemsCompressed && *GlobalBestCost == *OldGlobalBestCost)
            Alg->SolveCompressedSubproblem(*CurrentSubproblem, *Subproblems,
                                      &*GlobalBestCost);
    } else {
        int mid = (start + end) / 2;
        KarpPartition(start, mid,Alg);
        KarpPartition(mid + 1, end,Alg);
    }
}

/*
 * The CalculateSubproblems function is used to calculate the number of 
 * subproblems (Subproblems) created by the Karpartition function.
 */

static void CalculateSubproblems(int start, int end, LKH::LKHAlg * Alg)
{
    if (end - start + 1 <= Alg->SubproblemSize)
        (*Subproblems)++;
    else {
        int mid = (start + end) / 2;
        CalculateSubproblems(start, mid,Alg);
        CalculateSubproblems(mid + 1, end,Alg);
    }
}
