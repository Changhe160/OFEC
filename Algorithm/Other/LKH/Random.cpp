/*
 * This file contains a portable random generator. It will give
 * identical sequences of random integers for any platform with
 * at least 32-bit integers.
 *
 * A version of this generator is described in J. Bentley's column, 
 * "The Software Exploratorium", Unix Review 1991. It is based on 
 * Algorithm A in D. E. Knuth, The Art of Computer Programming, 
 * Vol 2, Section 3.2.2, pp. 172.  
 *  
 * The Random function returns a pseudo-random integer in the range
 * 0...INT_MAX-1.
 *   
 * The SRandom function uses the given seed for a new sequence of
 * pseudo-random numbers.  
 */

#include "LKH.h"
#include "../../../Global/global.h"

#undef STDLIB_RANDOM
/* #define STDLIB_RANDOM */

#ifdef STDLIB_RANDOM
#include <stdlib.h>
unsigned LKH::LKHAlg::Random()
{
    return rand();
}

void LKH::LKHAlg::SRandom(unsigned Seed)
{
    srand(Seed);
}

#else

#include <limits.h>
#include "../../../Utility/include.h"
#define PRANDMAX INT_MAX

static thread_local unique_ptr<int> a, b, initialized;
static thread_local unique_ptr<vector<int> > arr;

unsigned LKH::LKHAlg::Random()
{
	if(!a.get())
	{
		a.reset(new int(0));
		b.reset(new int(24));
		initialized.reset(new int(0));
		arr.reset(new vector<int>(55));
	}

	return Global::msp_global->getRandInt(0,PRANDMAX,Program_Problem);
  /*  int t;

    if (!*initialized)
        SRandom(7913);
    if ((*a)-- == 0)
        *a = 54;
    if ((*b)-- == 0)
        *b = 54;
	if ((t = arr.get()[*a] - arr.get()[*b]) < 0)
        t += PRANDMAX;
    return (arr.get()[*a] = t); */
}

void LKH::LKHAlg::SRandom(unsigned Seed)
{
    int i, ii, last, next;
	if(!a.get())
	{
		a.reset(new int(0));
		b.reset(new int(24));
		initialized.reset(new int(0));
		arr.reset(new vector<int>(55));
	}
    Seed %= PRANDMAX;
    (*arr.get())[0] = last = Seed;
    for (next = i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        (*arr.get())[ii] = next;
        if ((next = last - next) < 0)
            next += PRANDMAX;
        last = (*arr.get())[ii];
    }
    *initialized = 1;
    *a = 0;
    *b = 24;
    for (i = 0; i < 165; i++)
        Random();
}

#endif
