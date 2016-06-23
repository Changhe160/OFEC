#ifndef _GENETIC_H
#define _GENETIC_H

/*
 * This header specifies the interface for the genetic algorithm part of LKH.
 */
#include "../../../Utility/include.h"

typedef void (*CrossoverFunction) (LKH::LKHAlg *Alg);

extern thread_local unique_ptr<int> MaxPopulationSize; /* The maximum size of the population */ 
extern thread_local unique_ptr<int> PopulationSize;    /* The current size of the population */ 
extern thread_local unique_ptr<CrossoverFunction> Crossover;
extern thread_local unique_ptr<int*> Population;      /* Array of individuals (solution tours) */
extern thread_local unique_ptr<GainType> Fitness;     /* The fitness (tour cost) of each individual */

void AddToPopulation(GainType Cost,LKH::LKHAlg *Alg);
void ApplyCrossover(int i, int j,LKH::LKHAlg *Alg);
void FreePopulation();
int HasFitness(GainType Cost);
int LinearSelection(int Size, double Bias,LKH::LKHAlg *Alg);
GainType MergeTourWithIndividual(int i,LKH::LKHAlg *Alg);
void PrintPopulation(LKH::LKHAlg *Alg);
void ReplaceIndividualWithTour(int i, GainType Cost,LKH::LKHAlg *Alg);
int ReplacementIndividual(GainType Cost,LKH::LKHAlg *Alg);

void ERXT(LKH::LKHAlg *Alg);

#endif
