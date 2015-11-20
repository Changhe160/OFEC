#include "LKH.h"
#include "Delaunay.h"

/*
 * The SolveDelaunaySubproblems function attempts to improve a given tour
 * by means of a partitioning scheme based on Delaunay clustering. 
 *
 * The set of nodes is subdivided into clusters, with at most SubproblemSize 
 * nodes in each one. Each cluster together with the given tour induces a 
 * subproblem consisting of all nodes in the cluster, and with edges fixed 
 * between nodes that are connected by tour segments whose interior points
 * do not belong the cluster.  
 *  
 * If an improvement is found, the new tour is written to TourFile. 
 * The original tour is given by the SubproblemSuc references of the nodes.
 */

static int DelaunayClustering(int MaxClusterSize ,LKH::LKHAlg *Alg);

void LKH::LKHAlg::SolveDelaunaySubproblems()
{
    Node *N;
    GainType GlobalBestCost, OldGlobalBestCost;
    int CurrentSubproblem, Subproblems;
    double EntryTime = GetTime();

    AllocateStructures();
    ReadPenalties();

    /* Compute upper bound for the original problem */
    GlobalBestCost = 0;
    N = FirstNode;
    do {
        if (!Fixed(N, N->SubproblemSuc))
            GlobalBestCost += (this->*Distance)(N, N->SubproblemSuc);
        N->Subproblem = 0;
    }
    while ((N = N->SubproblemSuc) != FirstNode);
    if (TraceLevel >= 1) {
        if (TraceLevel >= 2)
            printff("\n");
        printff("*** Delaunay partitioning *** [Cost = " GainFormat "]\n",
                GlobalBestCost);
    }
    Subproblems = DelaunayClustering(SubproblemSize,this);
    for (CurrentSubproblem = 1;
         CurrentSubproblem <= Subproblems; CurrentSubproblem++) {
        OldGlobalBestCost = GlobalBestCost;
        SolveSubproblem(CurrentSubproblem, Subproblems, &GlobalBestCost);
        if (SubproblemsCompressed && GlobalBestCost == OldGlobalBestCost)
            SolveCompressedSubproblem(CurrentSubproblem, Subproblems,
                                      &GlobalBestCost);
    }
    printff("\nCost = " GainFormat, GlobalBestCost);
    if (Optimum != MINUS_INFINITY && Optimum != 0)
        printff(", Gap = %0.4f%%",
                100.0 * (GlobalBestCost - Optimum) / Optimum);
    printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
            GlobalBestCost < Optimum ? "<" : GlobalBestCost ==
            Optimum ? "=" : "");
    if (SubproblemBorders && Subproblems > 1)
        SolveSubproblemBorderProblems(Subproblems, &GlobalBestCost);
}

typedef struct Edge {
    LKH::LKHAlg::Node *From, *To;
    int Cost;
} Edge;

static int compareCost(const void *e1, const void *e2)
{
    int d1 = ((Edge *) e1)->Cost;
    int d2 = ((Edge *) e2)->Cost;
    return d1 < d2 ? -1 : d1 == d2 ? 0 : 1;
}

static int compareFromTo(const void *e1, const void *e2)
{
    int From1 = ((Edge *) e1)->From->Id;
    int To1 = ((Edge *) e1)->To->Id;
    int From2 = ((Edge *) e2)->From->Id;
    int To2 = ((Edge *) e2)->To->Id;
    return From1 < From2 ? -1 : From1 > From2 ? 1 : To1 - To2;
}

#define Size V

static int DelaunayClustering(int MaxClusterSize,LKH::LKHAlg *Alg)
{
    int Count = 0, Count1, Count2 = 0, i, j = 0;
    Edge *EdgeSet, Key;
    LKH::LKHAlg::Node *From, *To, *N, *F, *T;
    point *u, *v;
    edge *e_start, *e;

	delaunay(Alg->Dimension,Alg);
    for (i = 0; i < Alg->Dimension; i++) {
        u = &p_array.get()[i];
        e_start = e = u->entry_pt;
        do {
            v = Other_point(e, u);
            if (u->id < v->id)
                Count++;
        } while ((e = Next(e, u)) != e_start);
    }
    assert(EdgeSet = (Edge *) malloc(Count * sizeof(Edge)));
    for (i = 0; i < Alg->Dimension; i++) {
        u = &p_array.get()[i];
        e_start = e = u->entry_pt;
        do {
            v = Other_point(e, u);
            if (u->id < v->id) {
                EdgeSet[j].From =  From = &Alg->NodeSet[u->id];;
                EdgeSet[j].To = To = &Alg->NodeSet[v->id];;
                EdgeSet[j++].Cost = (Fixed(From, To) || Alg->IsCommonEdge(From, To)) ? INT_MIN :
                    (Alg->*(Alg->Distance))(From, To) * Alg->Precision + From->Pi + To->Pi;
            }
        } while ((e = Next(e, u)) != e_start);
    }
    free_memory();
    if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEOM ||
        Alg->WeightType == LKH::GEO_MEEUS || Alg->WeightType == LKH::GEOM_MEEUS) {
        N = Alg->FirstNode;
        while ((N = N->Suc) != Alg->FirstNode)
            if ((N->Y > 0) != (Alg->FirstNode->Y > 0))
                break;
        if (N != Alg->FirstNode) {
            N = Alg->FirstNode;
            do
                N->Zc = N->Y;
            while ((N = N->Suc) != Alg->FirstNode);
            /* Transform longitude (180 and -180 map to 0) */
            From = Alg->FirstNode;
            do {
                From->Zc = From->Y;
                if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
                    From->Y =
                        (int) From->Y + 5.0 * (From->Y -
                                               (int) From->Y) / 3.0;
                From->Y += From->Y > 0 ? -180 : 180;
                if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
                    From->Y =
                        (int) From->Y + 3.0 * (From->Y -
                                               (int) From->Y) / 5.0;
            } while ((From = From->Suc) != Alg->FirstNode);
			delaunay(Alg->Dimension,Alg);
            do
                From->Y = From->Zc;
            while ((From = From->Suc) != Alg->FirstNode);

            qsort(EdgeSet, Count, sizeof(Edge), compareFromTo);
            for (i = 0; i < Alg->Dimension; i++) {
                u = &p_array.get()[i];
                e_start = e = u->entry_pt;
                do {
                    v = Other_point(e, u);
                    if (u->id < v->id)
                        Count2++;
                } while ((e = Next(e, u)) != e_start);
            }
            Count1 = Count;
            assert(EdgeSet =
                   (Edge *) realloc(EdgeSet,
                                    (Count1 + Count2) * sizeof(Edge)));
            for (i = 0; i < Alg->Dimension; i++) {
                u = &p_array.get()[i];
                e_start = e = u->entry_pt;
                do {
                    v = Other_point(e, u);
                    if (u->id > v->id)
                        continue;
                    Key.From = From = &Alg->NodeSet[u->id];
                    Key.To = To = &Alg->NodeSet[v->id];
                    if (!bsearch
                        (&Key, EdgeSet, Count1, sizeof(Edge),
                         compareFromTo)) {
                        EdgeSet[Count].From = From;
                        EdgeSet[Count].To = To;
                        EdgeSet[Count].Cost =
                            (Fixed(From, To) || Alg->IsCommonEdge(From, To)) ? INT_MIN :
                           (Alg->*(Alg->Distance))(From,
                                     To) * Alg->Precision + From->Pi + To->Pi;
                        Count++;
                    }
                } while ((e = Next(e, u)) != e_start);
            }
            free_memory();
        }
    }
    qsort(EdgeSet, Count, sizeof(Edge), compareCost);

    /* Union-Find with path compression */
    N = Alg->FirstNode;
    do {
        N->Next = N;
        N->Size = 1;
    } while ((N = N->Suc) != Alg->FirstNode);
    for (i = 0; i < Count; i++) {
        for (F = EdgeSet[i].From; F != F->Next;)
            F = F->Next = F->Next->Next;
        for (T = EdgeSet[i].To; T != T->Next;)
            T = T->Next = T->Next->Next;
        if (F != T && F->Size + T->Size <= MaxClusterSize) {
            if (F->Size < T->Size) {
                F->Next = T;
                T->Size += F->Size;
            } else {
                T->Next = F;
                F->Size += T->Size;
            }
        }
    }
    free(EdgeSet);

    Count = 0;
    N = Alg->FirstNode;
    do {
        if (N->Next == N && N->Size > 3)
            N->Subproblem = ++Count;
    } while ((N = N->Suc) != Alg->FirstNode);
    do {
        for (T = N; T != T->Next;)
            T = T->Next = T->Next->Next;
        N->Subproblem = T->Subproblem;
    } while ((N = N->Suc) != Alg->FirstNode);
    return Count;
}
