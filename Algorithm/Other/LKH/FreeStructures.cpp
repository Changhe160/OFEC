#include "LKH.h"
#include "Sequence.h"
#include "Genetic.h"

/*      
 * The FreeStructures function frees all allocated structures.
 */

#define Free(s) { free(s); s = 0; }

void LKH::LKHAlg::FreeStructures()
{
    FreeCandidateSets();
    FreeSegments();
    if (NodeSet) {
        int i;
        for (i = 1; i <= Dimension; i++) {
            Node *N = &NodeSet[i];
            Free(N->MergeSuc);
            N->C = 0;
        }
        Free(NodeSet);
    }
    Free(CostMatrix);
    Free(BestTour);
    Free(BetterTour);
    Free(SwapStack);
    Free(HTable);
    Free(Rand);
    Free(CacheSig);
    Free(CacheVal);
    Free(Name);
    Free(Type);
    Free(EdgeWeightType);
    Free(EdgeWeightFormat);
    Free(EdgeDataFormat);
    Free(NodeCoordType);
    Free(DisplayDataType);
    Free(Heap);
    t.reset();
    T.reset();
    tSaved.reset();
    p.reset();
    q.reset();
    incl.reset();
    cycle.reset();
    G.reset();
    FreePopulation();
}

/*      
   The FreeSegments function frees the segments.
 */

void LKH::LKHAlg::FreeSegments()
{
    if (FirstSegment) {
        Segment *S = FirstSegment, *SPrev;
        do {
            SPrev = S->Pred;
            Free(S);
        }
        while ((S = SPrev) != FirstSegment);
        FirstSegment = 0;
    }
    if (FirstSSegment) {
        SSegment *SS = FirstSSegment, *SSPrev;
        do {
            SSPrev = SS->Pred;
            Free(SS);
        }
        while ((SS = SSPrev) != FirstSSegment);
        FirstSSegment = 0;
    }
}

/*      
 * The FreeCandidateSets function frees the candidate sets.
 */

void LKH::LKHAlg::FreeCandidateSets()
{
    Node *N = FirstNode;
    if (!N)
        return;
    do {
        Free(N->CandidateSet);
        Free(N->BackboneCandidateSet);
    }
    while ((N = N->Suc) != FirstNode);
}
