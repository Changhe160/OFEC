#ifndef _HEAP_H
#define _HEAP_H

/* 
 * This header specifies the interface for the use of heaps. 
 */

#include "LKH.h"

void MakeHeap(int Size,LKH::LKHAlg *Alg);
void HeapInsert(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);
void HeapDelete(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);
LKH::LKHAlg::Node *HeapDeleteMin(LKH::LKHAlg *Alg);
void HeapLazyInsert(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);
void Heapify(LKH::LKHAlg *Alg);
void HeapSiftUp(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);
void HeapSiftDown(LKH::LKHAlg::Node * N,LKH::LKHAlg *Alg);

#endif
