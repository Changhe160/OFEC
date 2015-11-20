#include "Segment.h"
#include "LKH.h"
#include "Heap.h"
#include "Sequence.h"

/*      
 * The AllocateStructures function allocates all necessary 
 * structures except nodes and candidates.
 */

#define Free(s) { free(s); s = 0; }

void LKH::LKHAlg::AllocateStructures()
{
	if(!K.get())
		K.reset(new int(0));
    int i, K;

    Free(Heap);
    Free(BestTour);
    Free(BetterTour);
    Free(HTable);
    Free(Rand);
    Free(CacheSig);
    Free(CacheVal);
	T.reset();
    G.reset();
    t.reset();
    p.reset();
    q.reset();
    Free(SwapStack);
    tSaved.reset();

    MakeHeap(Dimension,this);
    assert(BestTour = (int *) calloc(1 + Dimension, sizeof(int)));
    assert(BetterTour = (int *) calloc(1 + Dimension, sizeof(int)));
    assert(HTable = (HashTable *) malloc(sizeof(HashTable)));
    HashInitialize((HashTable *) HTable);
    SRandom(Seed);
    assert(Rand = (unsigned *)
           malloc((Dimension + 1) * sizeof(unsigned)));
    for (i = 1; i <= Dimension; i++)
        Rand[i] = Random();
    SRandom(Seed);
    if (WeightType != EXPLICIT) {
        for (i = 0; (1 << i) < (Dimension << 1); i++);
        i = 1 << i;
        assert(CacheSig = (int *) calloc(i, sizeof(int)));
        assert(CacheVal = (int *) calloc(i, sizeof(int)));
        CacheMask = i - 1;
    }
    AllocateSegments();
    K = MoveType;
    if (SubsequentMoveType > K)
        K = SubsequentMoveType;
	T.reset( new vector<Node *>(1 + 2 * K));
	G.reset(new vector<GainType>(2 * K));
	t.reset(new vector<Node *>(6 * K));
	tSaved.reset(new vector<Node *>(1 + 2 * K));
	p.reset(new vector<int>(6 * K));
	q.reset(new vector<int>(6 * K));
	incl.reset(new vector<int>(6 * K));
	cycle.reset(new vector<int>(6 * K));
    assert(SwapStack =
           (SwapRecord *) malloc((MaxSwaps + 6 * K) * sizeof(SwapRecord)));
}

/*      
 * The AllocateSegments function allocates the segments of the two-level tree.
 */

void LKH::LKHAlg::AllocateSegments()
{
    Segment *S = 0, *SPrev;
    SSegment *SS = 0, *SSPrev;
    int i;

    FreeSegments();
#ifdef THREE_LEVEL_TREE
    GroupSize = pow((double) Dimension, 1.0 / 3.0);
#elif defined TWO_LEVEL_TREE
    GroupSize = sqrt((double) Dimension);
#else
    GroupSize = Dimension;
#endif
    Groups = 0;
    for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S) {
        assert(S = (Segment *) malloc(sizeof(Segment)));
        S->Rank = ++Groups;
        if (!SPrev)
            FirstSegment = S;
        else
            SLink(SPrev, S);
    }
    SLink(S, FirstSegment);
#ifdef THREE_LEVEL_TREE
    SGroupSize = sqrt((double) Groups);
#else
    SGroupSize = Dimension;
#endif
    SGroups = 0;
    for (i = Groups, SSPrev = 0; i > 0; i -= SGroupSize, SSPrev = SS) {
        SS = (SSegment *) malloc(sizeof(SSegment));
        SS->Rank = ++SGroups;
        if (!SSPrev)
            FirstSSegment = SS;
        else
            SLink(SSPrev, SS);
    }
    SLink(SS, FirstSSegment);
}
