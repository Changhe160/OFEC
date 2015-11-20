#include "LKH.h"

/*
 * The KSwapKick function makes a random K-swap kick, K>=3 
 * (a generalization of the double-bridge kick).
 *
 * The algorithm is inspired by the thesis 
 *
 *    D. Richter,
 *    Toleranzen in Helsgauns Lin-Kernighan.Heuristik fur das TSP,
 *    Diplomarbeit, Martin-Luther-Universitat Halle-Wittenberg, 2006.
 */

static LKH::LKHAlg::Node *RandomNode(LKH::LKHAlg * Alg);
static int compare(const void *Na, const void *Nb);

void LKH::LKHAlg::KSwapKick(int K)
{
    Node **s, *N;
    int Count, i;

    assert(s = (Node **) malloc(K * sizeof(Node *)));
    Count = 0;
    N = FirstNode;
    do {
        N->Rank = ++Count;
        N->V = 0;
    } while ((N = N->Suc) != FirstNode);
    N = s[0] = RandomNode(this);
    if (!N)
        goto End_KSwapKick;
    N->V = 1;
    for (i = 1; i < K; i++) {
        N = s[i] = RandomNode(this);
        if (!N)
            K = i;
        else
            N->V = 1;
    }
    if (K < 3)
        goto End_KSwapKick;
    qsort(s, K, sizeof(Node *), compare);
    for (i = 0; i < K; i++)
        s[i]->OldSuc = s[i]->Suc;
    for (i = 0; i < K; i++)
        Link(s[(i + 2) % K], s[i]->OldSuc);
  End_KSwapKick:
    free(s);
}

/*
 * The RandomNode function returns a random node N, for
 * which the edge (N, N->Suc) is neither a fixed edge nor
 * a common edge of tours to be merged, and N has not 
 * previously been chosen.
 */

static LKH::LKHAlg::Node *RandomNode(LKH::LKHAlg * Alg)
{
    LKH::LKHAlg::Node *N;
    int Count;

    if (Alg->Dimension == Alg->DimensionSaved)
        N = &Alg->NodeSet[1 + Alg->Random() % Alg->Dimension];
    else {
        N = Alg->FirstNode;
        for (Count = Alg->Random() % Alg->Dimension; Count > 0; Count--)
            N = N->Suc;
    }
    Count = 0;
    while (((Fixed(N, N->Suc) || Alg->IsCommonEdge(N, N->Suc)) || N->V) && Count < Alg->Dimension) {
        N = N->Suc;
        Count++;
    }
    return Count < Alg->Dimension ? N : 0;
}

static int compare(const void *Na, const void *Nb)
{
    return (*(LKH::LKHAlg::Node **) Na)->Rank - (*(LKH::LKHAlg::Node **) Nb)->Rank;
}
