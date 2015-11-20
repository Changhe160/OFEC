#include "Segment.h"
#include "LKH.h"
#include "Sequence.h"

/*
 * The PatchCycles function attempts to improve the tour by patching the 
 * M >= 2 cycles that would appear if the move defined by t[1..2k] and 
 * incl[1..2k] was made. If the composite move results in a shorter 
 * tour, then the move is made, and the function returns the gain.
 *
 * On entry, Gain is the gain that could be obtained by making the non-
 * feasible move defined by t[1..2k] and incl[1..2k].
 *   
 * The function tries to patch the cycles by interleaving the alternating 
 * path represented by t with one or more alternating cycles. 
 *   
 * The function is called from BestKOptMove.   
 */

static GainType PatchCyclesRec(int k, int m, int M, GainType G0, LKH::LKHAlg *Alg);
static int ShortestCycle(int M, int k,LKH::LKHAlg *Alg);
static int Cycle(LKH::LKHAlg::Node * N, int k,LKH::LKHAlg *Alg);

static boost::thread_specific_ptr<int> CurrentCycle, Patchwork, RecLevel;

/*
 * The PatchCycles function tries to find a gainful move by patch the cycles 
 * that would occur if the move represented by t[1..2k] and incl[1..2k] was 
 * made using one one or more alternating cycles.
 * The alternating cycles are put in continuation of t, starting at 2k+1.
 */
void LKH::LKHAlg::freePatchCycle()
{
	CurrentCycle.reset();
	Patchwork.reset();
	RecLevel.reset();
}

GainType LKH::LKHAlg::PatchCycles(int k, GainType Gain)
{
    Node *s1, *s2, *sStart, *sStop;
    GainType NewGain;
    int M, i;

	if(!CurrentCycle.get())
	{
		CurrentCycle.reset(new int(0));
		Patchwork.reset(new int(0));
		RecLevel.reset(new int(0));
	}
	FindPermutation(k,this);
    M = Cycles(k);
    if (M == 1 && Gain > 0) {
        MakeKOptMove(k);
        return Gain;
    }
    if (M == 1 || M > PatchingC || k + M > NonsequentialMoveType)
        return 0;
    if (*RecLevel == 0)
        *Patchwork = 0;
    *CurrentCycle = ShortestCycle(M, k,this);
    for (i = 0; i < k; i++) {
        if ((*cycle.get())[(*p.get())[2 * i]] != *CurrentCycle)
            continue;
        sStart = (*t.get())[(*p.get())[2 * i]];
        sStop = (*t.get())[(*p.get())[2 * i + 1]];
        for (s1 = sStart; s1 != sStop; s1 = s2) {
            s2 = SUC(s1);
            if (FixedOrCommon(s1, s2))
                continue;
            if (++(*Patchwork) > Dimension)
                return 0;
            (*t.get())[2 * k + 1] = s1;
            (*t.get())[2 * k + 2] = s2;
            MarkDeleted(s1, s2);
            /* Find a set of gainful alternating cycles */
            NewGain = PatchCyclesRec(k, 2, M, Gain + (this->*C)(s1, s2),this);
            UnmarkDeleted(s1, s2);
            if (NewGain > 0)
                return NewGain;
        }
    }
    return 0;
}

static GainType PatchCyclesRec(int k, int m, int M, GainType G0, LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *s1, *s2, *s3, *s4, *s5, *s6, *S3 = 0, *S4 = 0;
    LKH::LKHAlg::Candidate *Ns2, *Ns4;
    GainType G1, G2, G3, G4, Gain, CloseUpGain,
        BestCloseUpGain = Alg->PatchingAExtended ? MINUS_INFINITY : 0;
    int X4, X6;
    int i, NewCycle, *cycleSaved = 0, *pSaved = 0;
    int Breadth2 = 0, Breadth4;

    s1 = (*t.get())[2 * k + 1];
    s2 = (*t.get())[i = 2 * (k + m) - 2];
    (*incl.get())[(*incl.get())[i] = i + 1] = i;

    /* Choose (s2,s3) as a candidate edge emanating from s2 */
    for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
        if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3) ||
			(NewCycle = Cycle(s3, k,Alg)) == *CurrentCycle)
            continue;
        if (++Breadth2 > Alg->MaxBreadth)
            break;
        MarkAdded(s2, s3);
        (*t.get())[2 * (k + m) - 1] = s3;
        G1 = G0 - Ns2->Cost;
        /* Choose s4 as one of s3's two neighbors on the tour */
        for (X4 = 1; X4 <= 2; X4++) {
            s4 = X4 == 1 ? s3->Pred : s3->Suc;
            if ((Fixed(s3, s4) || Alg->IsCommonEdge(s3, s4)) || Deleted(s3, s4))
                continue;
            MarkDeleted(s3, s4);
            (*t.get())[2 * (k + m)] = s4;
			G2 = G1 + (Alg->*(Alg->C))(s3, s4);
            if (M > 2) {
                if (!cycleSaved) {
                    assert(cycleSaved =
                           (int *) malloc(2 * k * sizeof(int)));
                    memcpy(cycleSaved, cycle.get() + 1, 2 * k * sizeof(int));
                }
                for (i = 1; i <= 2 * k; i++)
                    if ((*cycle.get())[i] == NewCycle)
                        (*cycle.get())[i] = *CurrentCycle;
                /* Extend the current alternating path */
				if ((Gain = PatchCyclesRec(k, m + 1, M - 1, G2,Alg)) > 0) {
                    UnmarkAdded(s2, s3);
                    UnmarkDeleted(s3, s4);
                    goto End_PatchCyclesRec;
                }
                memcpy(cycle.get() + 1, cycleSaved, 2 * k * sizeof(int));
                if (Alg->PatchingA >= 2 && *Patchwork < Alg->Dimension &&
                    k + M < Alg->NonsequentialMoveType &&
                    !Alg->Forbidden(s4, s1) &&
                    (!Alg->PatchingARestricted || Alg->IsCandidate(s4, s1))) {
                    GainType Bound = BestCloseUpGain >= 0 ||
                        Alg->IsCandidate(s4, s1) ? BestCloseUpGain : 0;
                    if ((!Alg->c || G2 - (Alg->*(Alg->c))(s4, s1) > Bound) &&
                        (CloseUpGain = G2 - (Alg->*(Alg->C))(s4, s1)) > Bound) {
                        S3 = s3;
                        S4 = s4;
                        BestCloseUpGain = CloseUpGain;
                    }
                }
            } else if (!Alg->Forbidden(s4, s1) && (!Alg->c || G2 - (Alg->*(Alg->c))(s4, s1) > 0)
                       && (Gain = G2 - (Alg->*(Alg->C))(s4, s1)) > 0) {
                (*incl.get())[(*incl.get())[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
                Alg->MakeKOptMove(k + m);
                UnmarkAdded(s2, s3);
                UnmarkDeleted(s3, s4);
                goto End_PatchCyclesRec;
            }
            UnmarkDeleted(s3, s4);
        }
        UnmarkAdded(s2, s3);
    }
    if (M == 2 && !Alg->PatchingCRestricted) {
        /* Try to patch the two cycles by a sequential 3-opt move */
        (*incl.get())[(*incl.get())[2 * (k + m)] = 2 * (k + m) + 1] = 2 * (k + m);
        (*incl.get())[(*incl.get())[2 * k + 1] = 2 * (k + m) + 2] = 2 * k + 1;
        Breadth2 = 0;
        /* Choose (s2,s3) as a candidate edge emanating from s2 */
        for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
            if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3))
                continue;
            if (++Breadth2 > Alg->MaxBreadth)
                break;
            (*t.get())[2 * (k + m) - 1] = s3;
            G1 = G0 - Ns2->Cost;
			NewCycle = Cycle(s3, k,Alg);
            /* Choose s4 as one of s3's two neighbors on the tour */
            for (X4 = 1; X4 <= 2; X4++) {
                s4 = X4 == 1 ? s3->Pred : s3->Suc;
                if ((Fixed(s3, s4) || Alg->IsCommonEdge(s3, s4)) || Deleted(s3, s4))
                    continue;
                (*t.get())[2 * (k + m)] = s4;
                G2 = G1 + (Alg->*(Alg->C))(s3, s4);
                Breadth4 = 0;
                /* Choose (s4,s5) as a candidate edge emanating from s4 */
                for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
                    if (s5 == s4->Pred || s5 == s4->Suc || s5 == s1 ||
                        Added(s4, s5) ||
                        (NewCycle == *CurrentCycle &&
						Cycle(s5, k,Alg) == *CurrentCycle))
                        continue;
                    if (++Breadth4 > Alg->MaxBreadth)
                        break;
                    G3 = G2 - Ns4->Cost;
                    /* Choose s6 as one of s5's two neighbors on the tour */
                    for (X6 = 1; X6 <= 2; X6++) {
                        s6 = X6 == 1 ? s5->Pred : s5->Suc;
                        if (s6 == s1 || Alg->Forbidden(s6, s1)
                            || (Fixed(s5, s6) || Alg->IsCommonEdge(s5, s6))
                            || Deleted(s5, s6)
                            || Added(s6, s1))
                            continue;
                        G4 = G3 + (Alg->*(Alg->C))(s5, s6);
                        if ((!Alg->c || G4 - (Alg->*(Alg->c))(s6, s1) > 0) &&
                            (Gain = G4 - (Alg->*(Alg->C))(s6, s1)) > 0) {
                            if (!pSaved) {
                                assert(pSaved =
                                       (int *) malloc(2 * k *
                                                      sizeof(int)));
                                memcpy(pSaved, p.get() + 1, 2 * k * sizeof(int));
                            }
                            (*t.get())[2 * (k + m) + 1] = s5;
                            (*t.get())[2 * (k + m) + 2] = s6;
                            if (FeasibleKOptMove(k + m + 1)) {
                                Alg->MakeKOptMove(k + m + 1);
                                goto End_PatchCyclesRec;
                            }
                            memcpy(p.get() + 1, pSaved, 2 * k * sizeof(int));
                            for (i = 1; i <= 2 * k; i++)
                                (*q.get())[(*p.get())[i]] = i;
                        }
                    }
                }
            }
        }
    }
    Gain = 0;
    if (S4) {
        int OldCycle = *CurrentCycle;
        if (!pSaved) {
            assert(pSaved = (int *) malloc(2 * k * sizeof(int)));
            memcpy(pSaved, p.get() + 1, 2 * k * sizeof(int));
        }
        (*t.get())[2 * (k + m) - 1] = S3;
        (*t.get())[2 * (k + m)] = S4;
        (*incl.get())[(*incl.get())[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
        /* Find a new alternating cycle */
        Alg->PatchingA--;
        (*RecLevel)++;
        MarkAdded(s2, S3);
        MarkDeleted(S3, S4);
        MarkAdded(S4, s1);
        Gain = Alg->PatchCycles(k + m, BestCloseUpGain);
        UnmarkAdded(s2, S3);
        UnmarkDeleted(S3, S4);
        UnmarkAdded(S4, s1);
        (*RecLevel)--;
        Alg->PatchingA++;
        if (Gain <= 0) {
            memcpy(cycle.get() + 1, cycleSaved, 2 * k * sizeof(int));
            memcpy(p.get() + 1, pSaved, 2 * k * sizeof(int));
            for (i = 1; i <= 2 * k; i++)
                (*q.get())[(*p.get())[i]] = i;
            *CurrentCycle = OldCycle;
        }
    }

  End_PatchCyclesRec:
    free(cycleSaved);
    free(pSaved);
    return Gain;
}

/*
 * The Cycle function returns the number of the cycle containing
 * a given node, N.
 *
 * Time complexity: O(log k). 
 */

static int Cycle(LKH::LKHAlg::Node * N, int k,LKH::LKHAlg *Alg)
{
    /* Binary search */
    int Low = 1, High = k;
    while (Low < High) {
        int Mid = (Low + High) / 2;
        if (Alg->Between_SL((*t.get())[(*p.get())[2 * Low]], N, (*t.get())[(*p.get())[2 * Mid + 1]]))
            High = Mid;
        else
            Low = Mid + 1;
    }
    return (*cycle.get())[(*p.get())[2 * Low]];
}

/*
 * The ShortestCycle function returns the number of the cycle with 
 * the smallest number of nodes. Note however that if the two-level 
 * list is used, the number of nodes of each cycle is only approximate
 * (for efficiency reasons). 
 *
 * Time complexity: O(k + M), where M = Cycles(k). 
 * 
 * The function may only be called after a call of the Cycles function.
 */

static int ShortestCycle(int M, int k,LKH::LKHAlg *Alg)
{
    int i, Cycle, MinCycle = 0;
    int *Size, MinSize = INT_MAX;

    assert(Size = (int *) calloc(1 + M, sizeof(int)));
    (*p.get())[0] = (*p.get())[2 * k];
    for (i = 0; i < 2 * k; i += 2)
        Size[(*cycle.get())[(*p.get())[i]]] += Alg->SegmentSize((*t.get())[(*p.get())[i]], (*t.get())[(*p.get())[i + 1]]);
    for (Cycle = 1; Cycle <= M; Cycle++) {
        if (Size[Cycle] < MinSize) {
            MinSize = Size[Cycle];
            MinCycle = Cycle;
        }
    }
    free(Size);
    return MinCycle;
}
