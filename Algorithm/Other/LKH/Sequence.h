#ifndef SEQUENCE_H
#define SEQUENCE_H

/* 
 * This header specifies the interface for the use of node sequences.
 *   
 * The functions BestKOptMove and BacktrackKOptMove are implemented 
 * by means of such sequences. 
 */

#include "LKH.h"

extern boost::thread_specific_ptr<vector<LKH::LKHAlg::Node *> > t;       /* The sequence of nodes to be used in a move */
extern boost::thread_specific_ptr<vector<LKH::LKHAlg::Node *> > T;       /* The sequence of nodes to be used in a move */
extern boost::thread_specific_ptr<vector<LKH::LKHAlg::Node *> > tSaved;       /* The sequence of nodes to be used in a move */
extern boost::thread_specific_ptr<vector<int> > p;         /* The permutation corresponding to the sequence in which the t's occur on the tour */
extern boost::thread_specific_ptr<vector<int> > q;         /* The inverse permutation of p */
extern boost::thread_specific_ptr<vector<int> > incl;      /* Array: incl[i] == j, if (t[i], t[j]) is an inclusion edge */
extern boost::thread_specific_ptr<vector<int> > cycle;     /* Array: cycle[i] is cycle number of t[i] */
extern boost::thread_specific_ptr<vector<GainType> > G;    /* For storing the G-values in the BestKOptMove function */
extern boost::thread_specific_ptr<int> K;          /* The value K for the current K-opt move */

int FeasibleKOptMove(int k);
void FindPermutation(int k,LKH::LKHAlg *Alg);
int Cycles(int k);

int Added(const LKH::LKHAlg::Node * ta, const LKH::LKHAlg::Node * tb);
int Deleted(const LKH::LKHAlg::Node * ta, const LKH::LKHAlg::Node * tb);

void MarkAdded(LKH::LKHAlg::Node *ta, LKH::LKHAlg::Node *tb);
void MarkDeleted(LKH::LKHAlg::Node *ta, LKH::LKHAlg::Node *tb);
void UnmarkAdded(LKH::LKHAlg::Node *ta, LKH::LKHAlg::Node *tb);
void UnmarkDeleted(LKH::LKHAlg::Node *ta, LKH::LKHAlg::Node *tb);

#endif
