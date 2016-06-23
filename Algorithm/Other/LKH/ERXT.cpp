#include "LKH.h"
#include "Genetic.h"

/* 
 * The ERXT function applies the Edge Recombination Crossover operator (ERX)
 * on the two individuals (tours) represented by the Suc and Next references,
 * resepectively.
 * 
 * ERX was originally described in 
 *
 *     D. Whitley, T. Starkweather, and D. Fuquay,
 *     Scheduling Problems and the Traveling Salesman:
 *     the Genetic Edge Recombination Operator. 
 *     Proc. Third Int. Conf. on Genetic Algorithms and Their Applications
 *     (1989)
 * 
 * ERXT implements the variant of ERX based on tabu-edges (Edge-T) described in
 *
 *     Chuan-Kang Ting,
 *     Improving Edge Recombination through Alternate Inheritance and
 *     Greedy Manner.
 *     Lecture Notes in Computer Science 3004, pp. 207-216 (2004)
 *
 * However, ERXT does not implement the greedy strategy used by Edge-T for 
 * choosing foreign edges.
 */
static thread_local unique_ptr<LKH::LKHAlg::Node> FirstFree;
static thread_local unique_ptr<int> Tabu;

static LKH::LKHAlg::Node *SelectNext(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);

void LKH::LKHAlg::freeERXT()
{
	FirstFree.release();
}

void ERXT(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *N, *Next;
    int i;
	if(!Tabu.get())
		Tabu.reset(new int(0));
    *Tabu = 0;
    N = Alg->FirstNode;
    do {
        N->OldSuc = N->Suc;
        N->OldSuc->OldPred = N;
        N->Next->Prev = N;
        N->Suc->Pred = N;
        N->V = 0;
    }
    while ((N = N->Suc) != Alg->FirstNode);
    if (Alg->Dimension == Alg->DimensionSaved)
        Alg->FirstNode = &Alg->NodeSet[1 + Alg->Random() % Alg->Dimension];
    else
        for (i = Alg->Random() % Alg->Dimension; i > 0; i--)
            Alg->FirstNode = Alg->FirstNode->Suc;
    N = Alg->FirstNode;
    N->V = 1;
	FirstFree.release();
	FirstFree.reset(N->Suc);
    N->Pred->Suc = N->Suc;
    N->Suc->Pred = N->Pred;
    for (i = 1; i < Alg->Dimension; i++) {
		Next = SelectNext(N,Alg);
        if (Next == FirstFree.get())
			FirstFree.reset(Next->Suc);
        Next->Pred->Suc = Next->Suc;
        Next->Suc->Pred = Next->Pred;
        Link(N, Next);
        N = Next;
        N->V = 1;
    }
    Link(N, Alg->FirstNode);
}

/*
 * The EdgeCount function computes the number of unused edges emanating
 * from a given node, N.
 */

static int EdgeCount(LKH::LKHAlg::Node * N)
{
    int Count = 0;
    LKH::LKHAlg::Node *Next;

    if (!N->OldPred->V)
        Count++;
    if (!N->OldSuc->V)
        Count++;
    Next = N->Prev;
    if (!Next->V && Next != N->OldPred && Next != N->OldSuc)
        Count++;
    Next = N->Next;
    if (!Next->V && Next != N->OldPred && Next != N->OldSuc)
        Count++;
    return Count;
}

#define IsCommonEdge(Na, Nb)\
    (((Na)->OldPred == (Nb) || (Na)->OldSuc == (Nb)) &&\
     ((Na)->Prev == (Nb) || (Na)->Next == (Nb)))

/*
 * The SelectNext function select the next node to be added as a neighbor
 * to a given node, N.
 *
 * The function chooses the neighbor node with the highest priority (See Ting's
 * paper). If two or more possible neighbor nodes have the same priority,
 * then one of them is chosen randomly. If the node has no neighbors on the two
 * two tours, then the first node in the list of unused nodes is chosen.
 */

static LKH::LKHAlg::Node *SelectNext(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Next, *Alternative[4];
    int Alternatives = 0, Score, MaxScore = INT_MIN, i;

    for (i = 1; i <= 4; i++) {
        Next = i == 1 ? N->OldPred : i == 2 ? N->OldSuc :
            i == 3 ? N->Prev : N->Next;
        if (!Next->V &&
            (i <= 2 || (Next != N->OldPred && Next != N->OldSuc))) {
            if (Fixed(N, Next))
                Score = INT_MAX;
            else {
                Score = IsCommonEdge(N, Next) ? 4 : 0;
                Score -= EdgeCount(Next);
                Score -= i <= 2 ? *Tabu : -*Tabu;
            }
            if (Score >= MaxScore) {
                if (Score > MaxScore)
                    Alternatives = 0;
                Alternative[Alternatives++] = Next;
                MaxScore = Score;
            }
        }
    }
    if (Alternatives > 0) {
        Next = Alternative[Alg->Random() % Alternatives];
        if (Next == N->OldPred || Next == N->OldSuc)
            (*Tabu)++;
        if (Next == N->Prev || Next == N->Next)
            (*Tabu)--;
        return Next;
    }
    return FirstFree.get();
}
