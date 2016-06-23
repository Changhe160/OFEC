#include "LKH.h"

static void BuildSubKDTree(int start, int end,LKH::LKHAlg * Alg);
static void Partition(int start, int end, int k, int axis);
static char FindMaxSpread(int start, int end, LKH::LKHAlg * Alg);
static void Swap(int i, int j);

static thread_local unique_ptr<LKH::LKHAlg::Node *> KDTree;
static thread_local unique_ptr<int> cutoff;

#define Coord(N, axis) (axis == 0 ? (N)->X : axis == 1 ? (N)->Y : (N)->Z)

/*
 * The BuildKDTree function builds a balanced K-d tree of all nodes.
 *
 * See
 *    Jon Louis Bentley: K-d Trees for Semidynamic Point Sets. 
 *    Symposium on Computational Geometry 1990: 187-197
 */

LKH::LKHAlg::Node ** LKH::LKHAlg::BuildKDTree(int Cutoff)
{
    int i;
    Node *N;
	if(!cutoff.get())
		cutoff.reset(new int(0));
    *cutoff = Cutoff >= 1 ? Cutoff : 1;
	if(!KDTree.get())
		KDTree.reset(new Node*[Dimension]);
    for (i = 0, N = FirstNode; i < Dimension; i++, N = N->Suc)
        KDTree.get()[i] = N;
    BuildSubKDTree(0, Dimension - 1,this);
    return KDTree.get();
}

/*
 * The BuildSubKDTree function arranges the nodes KDTree[start:end]
 * to leave a balanced K-d tree in KDTree[start:end]. 
 */

static void BuildSubKDTree(int start, int end, LKH::LKHAlg * Alg)
{
    if (end - start + 1 > *cutoff) {
        int mid = (start + end) / 2;
        char axis = FindMaxSpread(start, end, Alg);
        Partition(start, end, mid, axis);
        KDTree.get()[mid]->Axis = axis;
        BuildSubKDTree(start, mid - 1,Alg);
        BuildSubKDTree(mid + 1, end,Alg);
    }
}

/*
 * The Partition function partitions the K-d tree about the (k-1)th 
 * smallest element (the one in KDTree[k]): It arranges 
 * KDtree[start:end] to leave Coord(KDTree[start:k-1], axis) less 
 * than or equal to Coord(KDTree[k+1:end], axis).
 *
 * For example, we could call Partition(a, 0, N - 1, N/2, axis) to
 * partition KDTree on the median, leaving the median in KDTree[N/2].
 *
 * Partition is linear time on the avarage.
 */

static void Partition(int start, int end, int k, int axis)
{
    while (start < end) {
        int i = start, j = end - 1, mid = (start + end) / 2;
        double pivot;
        if (Coord(KDTree.get()[mid], axis) < Coord(KDTree.get()[start], axis))
            Swap(start, mid);
        if (Coord(KDTree.get()[end], axis) < Coord(KDTree.get()[start], axis))
            Swap(start, end);
        if (Coord(KDTree.get()[end], axis) < Coord(KDTree.get()[mid], axis))
            Swap(mid, end);
        if (end - start <= 2)
            return;
        Swap(mid, j);
        pivot = Coord(KDTree.get()[j], axis);
        while (1) {
            while (Coord(KDTree.get()[++i], axis) < pivot);
            while (pivot < Coord(KDTree.get()[--j], axis));
            if (i >= j)
                break;
            Swap(i, j);
        }
        Swap(i, end - 1);
        if (i >= k)
            end = i - 1;
        if (i <= k)
            start = i + 1;
    }
}

static void Swap(int i, int j)
{
    LKH::LKHAlg::Node *T = KDTree.get()[i];
    KDTree.get()[i] = KDTree.get()[j];
    KDTree.get()[j] = T;
}

/*
 * The FindMaxSpread returns the dimension with largest difference
 * between minimmum and maximum among the points in KDTree[start:end].
 */

static char FindMaxSpread(int start, int end,LKH::LKHAlg * Alg)
{
    int i, axis;
    LKH::LKHAlg::Node *N;
    double Min[3], Max[3];

    N = KDTree.get()[start];
    Min[0] = Max[0] = N->X;
    Min[1] = Max[1] = N->Y;
    Min[2] = Max[2] = N->Z;
    for (i = start + 1; i <= end; i++) {
        for (axis = Alg->CoordType == LKH::THREED_COORDS ? 2 : 1; axis >= 0; axis--) {
            N = KDTree.get()[i];
            if (Coord(N, axis) < Min[axis])
                Min[axis] = Coord(N, axis);
            else if (Coord(N, axis) > Max[axis])
                Max[axis] = Coord(N, axis);
        }
    }
    if (Max[0] - Min[0] > Max[1] - Min[1])
        return Alg->CoordType != LKH::THREED_COORDS
            || Max[0] - Min[0] > Max[2] - Min[2] ? 0 : 2;
    return Alg->CoordType != LKH::THREED_COORDS
        || Max[1] - Min[1] > Max[2] - Min[2] ? 1 : 2;
}
