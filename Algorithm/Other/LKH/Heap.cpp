#include "LKH.h"
#include "Heap.h"

/*
 * A binary heap is used to implement a priority queue. 
 *
 * A heap is useful in order to speed up the computations of minimum 
 * spanning trees. The elements of the heap are the nodes, and the
 * priorities (ranks) are their associated costs (their minimum distance 
 * to the current tree). 
 */
static boost::thread_specific_ptr<int> HeapCount,HeapCapacity;

/*      
 * The MakeHeap function creates an empty heap. 
 */

void MakeHeap(int Size,LKH::LKHAlg *Alg)
{
    assert(Alg->Heap = (LKH::LKHAlg::Node **) malloc((Size + 1) * sizeof(LKH::LKHAlg::Node *)));
	if(!HeapCount.get())
	{
		HeapCount.reset(new int(0));
		HeapCapacity.reset(new int(0));
	}
    *HeapCapacity = Size;
    *HeapCount = 0;
}

/*
 * The HeapSiftUp function is called when the rank of a node is decreased, 
 * or when a node is inserted into the heap.
 * The function moves the node forward in the heap (the foremost node
 * of the heap has the lowest rank).
 * When calling HeapSiftUp(N), node N must belong to the heap.              
 */

void HeapSiftUp(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
    int Loc = N->Loc, Parent = Loc / 2;

    while (Parent && N->Rank < Alg->Heap[Parent]->Rank) {
        Alg->Heap[Loc] = Alg->Heap[Parent];
        Alg->Heap[Loc]->Loc = Loc;
        Loc = Parent;
        Parent /= 2;
    }
    Alg->Heap[Loc] = N;
    N->Loc = Loc;
}

/*
 * The HeapSiftDown function is called by the Heapify and HeapDeleteMin 
 * functions. The function moves the node backwards in the heap 
 * (the foremost node of the heap has the lowest rank).
 * When calling HeapSiftDown(N), node N must belong to the heap.              
 */

void HeapSiftDown(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
    int Loc = N->Loc, Child;

    while (Loc <= *HeapCount / 2) {
        Child = 2 * Loc;
        if (Child < *HeapCount && Alg->Heap[Child + 1]->Rank < Alg->Heap[Child]->Rank)
            Child++;
        if (N->Rank <= Alg->Heap[Child]->Rank)
            break;
        Alg->Heap[Loc] = Alg->Heap[Child];
        Alg->Heap[Loc]->Loc = Loc;
        Loc = Child;
    }
    Alg->Heap[Loc] = N;
    N->Loc = Loc;
}

/*       
 * The HeapDeleteMin function deletes the foremost node from the heap. 
 * The function returns a pointer to the deleted node (0, if the heap
 * is empty).
 */

LKH::LKHAlg::Node *HeapDeleteMin(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Remove;

    if (!*HeapCount)
        return 0;
    Remove = Alg->Heap[1];
    Alg->Heap[1] = Alg->Heap[(*HeapCount)--];
    Alg->Heap[1]->Loc = 1;
	HeapSiftDown(Alg->Heap[1],Alg);
    Remove->Loc = 0;
    return Remove;
}

/*       
 * The HeapInsert function inserts a node N into the heap.
 * When calling HeapInsert(N), node N must not belong to the heap.
 */

void HeapInsert(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
	HeapLazyInsert(N,Alg);
	HeapSiftUp(N,Alg);
}

/*
 * The HeapDelete function deletes a node N from the heap.
 */

void HeapDelete(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
    int Loc = N->Loc;
    if (!Loc)
        return;
    Alg->Heap[Loc] = Alg->Heap[(*HeapCount)--];
    Alg->Heap[Loc]->Loc = Loc;
    if (Alg->Heap[Loc]->Rank > N->Rank)
		HeapSiftDown(Alg->Heap[Loc],Alg);
    else
		HeapSiftUp(Alg->Heap[Loc],Alg);
    N->Loc = 0;
}

/*       
 * The HeapLazyInsert function inserts a node as the last node of the heap.
 * This may destroy the heap condition, but it can later be restored by 
 * calling the Heapify function.
 * When calling HeapLazyInsert(N), node N must not belong to the heap.
 */

void HeapLazyInsert(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg)
{
    assert(*HeapCount < *HeapCapacity);
    Alg->Heap[++(*HeapCount)] = N;
    N->Loc = *HeapCount;
}

/*       
 * The Heapify function constructs a heap from its nodes.
 */

void Heapify(LKH::LKHAlg *Alg)
{
    int Loc;
    for (Loc = *HeapCount / 2; Loc >= 1; Loc--)
		HeapSiftDown(Alg->Heap[Loc],Alg);
}
