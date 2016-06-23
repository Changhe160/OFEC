#include "LKH.h"
#include "Genetic.h"

/*
 * The AddToPopulation function adds the current tour as an individual to 
 * the population. The fitness of the individual is set equal to the cost
 * of the tour. The population is kept sorted in increasing fitness order.
 */

thread_local unique_ptr<int> MaxPopulationSize; /* The maximum size of the population */ 
thread_local unique_ptr<int> PopulationSize;    /* The current size of the population */ 
thread_local unique_ptr<CrossoverFunction> Crossover;
thread_local unique_ptr<int*> Population;      /* Array of individuals (solution tours) */
thread_local unique_ptr<GainType> Fitness;     /* The fitness (tour cost) of each individual */

void AddToPopulation(GainType Cost,LKH::LKHAlg *Alg)
{
    int i, *P;
    LKH::LKHAlg::Node *N;
    if (!Population.get()) {
		Population.reset(new int*[*MaxPopulationSize]);
        for (i = 0; i < *MaxPopulationSize; i++)
            assert(Population.get()[i] =
                   (int *) malloc((1 + Alg->Dimension) * sizeof(int)));
		Fitness.reset(new GainType[*MaxPopulationSize]);
    }
    for (i = *PopulationSize; i >= 1 && Cost < Fitness.get()[i - 1]; i--) {
        Fitness.get()[i] = Fitness.get()[i - 1];
        P = Population.get()[i];
        Population.get()[i] = Population.get()[i - 1];
        Population.get()[i - 1] = P;
    }
    Fitness.get()[i] = Cost;
    P = Population.get()[i];
    N = Alg->FirstNode;
    i = 1;
    do
        P[i++] = N->Id;
    while ((N = N->Suc) != Alg->FirstNode);
    P[0] = P[Alg->Dimension];
    (*PopulationSize)++;
}

/*
 * The ApplyCrossover function applies a specified crossover operator to two 
 * individuals.
 */

void ApplyCrossover(int i, int j,LKH::LKHAlg *Alg)
{
    int *Pi, *Pj, k;

    Pi = Population.get()[i];
    Pj = Population.get()[j];
    for (k = 1; k <= Alg->Dimension; k++) {
        Alg->NodeSet[Pi[k - 1]].Suc = &Alg->NodeSet[Pi[k]];
        Alg->NodeSet[Pj[k - 1]].Next = &Alg->NodeSet[Pj[k]];
    }
    if (Alg->TraceLevel >= 1)
        Alg->printff("Crossover(%d,%d)\n", i + 1, j + 1);
    /* Apply the crossover operator */
	(*Crossover.get())(Alg);
}

#define Free(s) { free(s); s = 0; }

/*
 * The FreePopulation function frees the memory space allocated to the 
 * population.
 */

void FreePopulation()
{
    if (Population.get()) {
        int i;
        for (i = 0; i < *MaxPopulationSize; i++)
            Free(Population.get()[i]);
		Population.reset();
		Fitness.reset();
    }
    *PopulationSize = 0;
}

/*
 * The HasFitness function returns 1 if the population contains an
 * individual with fitness equal to a given tour cost; otherwise 0.
 *
 * Since the population is sorted in fitness order the test may be
 * made by binary search.
 */

int HasFitness(GainType Cost)
{
    int Low = 0, High = *PopulationSize - 1;
    while (Low < High) {
        int Mid = (Low + High) / 2;
        if (Fitness.get()[Mid] < Cost)
            Low = Mid + 1;
        else
            High = Mid;
    }
    return High >= 0 && Fitness.get()[High] == Cost;
}

/*
 * Random01 is an auxiliary function for computing a random double number
 * in the range [0;1).
 */

static double Random01(LKH::LKHAlg *Alg)
{
    return ((double) Alg->Random()) / INT_MAX;
}

/*
 * The LinearSelection function is used to select an individual with 
 * random linear bias towards the best members of the population.
 * The parameter Bias is a number between 1.0 and 2.0.
 *
 * See
 *     Darrell Whitley,
 *     The GENITOR algorithm and selection pressure:
 *     Why rank-based allocation of reproductive trials is best. 
 *     Proceedings of the Third International Conference on Genetic Algorithms,
 *     1989.
 */

int LinearSelection(int Size, double Bias,LKH::LKHAlg *Alg)
{
    return (int) (Size *
                  (Bias -
				  sqrt((Bias * Bias - 4 * (Bias - 1) * Random01(Alg)))) /
                  2 / (Bias - 1));
}

/*
 * The MergeTourWithIndividual function attempts to find a short tour by
 * merging the current tour with a specified inddividual of the population.
 * The merging algorithm is the iterative partial transcription algrithm
 * described by Mobius, Freisleben, Merz and Schreiber.
 */

GainType MergeTourWithIndividual(int i,LKH::LKHAlg *Alg)
{
    int *Pi, k;

    assert(i >= 0 && i < *PopulationSize);
    Pi = Population.get()[i];
    for (k = 1; k <= Alg->Dimension; k++)
        Alg->NodeSet[Pi[k - 1]].Next = &Alg->NodeSet[Pi[k]];
    return Alg->MergeWithTour();
}

/*
 * The PrintPopulation function prints the cost and gap to optimum for
 * each individual of the population.
 */

void PrintPopulation(LKH::LKHAlg *Alg)
{
    int i;
    Alg->printff("Population:\n");
    for (i = 0; i < *PopulationSize; i++) {
        Alg->printff("%3d: " GainFormat, i + 1, Fitness.get()[i]);
        if (Alg->Optimum != MINUS_INFINITY && Alg->Optimum != 0)
            Alg->printff(", Gap = %0.4f%%",
                    100.0 * (Fitness.get()[i] - Alg->Optimum) / Alg->Optimum);
        Alg->printff("\n");
    }
}

/*
 * The ReplaceIndividualWithTour function replaces a given individual in 
 * the population by an indidual that represents the current tour.
 * The population is kept sorted in increasing fitness order.
 */

void ReplaceIndividualWithTour(int i, GainType Cost,LKH::LKHAlg *Alg)
{
    int j, *P;
    LKH::LKHAlg::Node *N;

    assert(i >= 0 && i < *PopulationSize);
    Fitness.get()[i] = Cost;
    P = Population.get()[i];
    N = Alg->FirstNode;
    for (j = 1; j <= Alg->Dimension; j++) {
        P[j] = N->Id;
        N = N->Suc;
    }
    P[0] = P[Alg->Dimension];
    while (i >= 1 && Cost < Fitness.get()[i - 1]) {
        Fitness.get()[i] = Fitness.get()[i - 1];
        Population.get()[i] = Population.get()[i - 1];
        i--;
    }
    Fitness.get()[i] = Cost;
    Population.get()[i] = P;
    while (i < *PopulationSize - 1 && Cost > Fitness.get()[i + 1]) {
        Fitness.get()[i] = Fitness.get()[i + 1];
        Population.get()[i] = Population.get()[i + 1];
        i++;
    }
    Fitness.get()[i] = Cost;
    Population.get()[i] = P;
}

/* 
 * The DistanceToIndividual returns the number of different edges between 
 * the tour (given by OldSuc) and individual i. 
 */

static int DistanceToIndividual(int i,LKH::LKHAlg *Alg) { 
    int Count = 0, j, *P = Population.get()[i];
    LKH::LKHAlg::Node *N;
    
    for (j = 0; j < Alg->Dimension; j++) {
        N = &Alg->NodeSet[P[j]];
        (N->Next = &Alg->NodeSet[P[j + 1]])->Prev = N;
    }
    N = Alg->FirstNode;
    do
        if (N->OldSuc != N->Next && N->OldSuc != N->Prev)
            Count++;
    while ((N = N->OldSuc) != Alg->FirstNode);
    return Count;
}

/*
 * The ReplacementIndividual function returns the individual to be 
 * replaced with the current tour. The function implements the 
 * replacement strategy (CD/RW) proposed in
 *
 *      M. Lozano, F. Herrera, and J. R. Cano,
 *      Replacement strategies to preserve useful diversity in
 *      steady-state genetic algorithms.
 *      Information Sciences 178 (2008) 4421–4433.
 */

int ReplacementIndividual(GainType Cost,LKH::LKHAlg *Alg) {
    int i, j, d, *P;
    int MinDist = INT_MAX, CMin = *PopulationSize - 1;
    LKH::LKHAlg::Node *N = Alg->FirstNode;
    while ((N = N->OldSuc = N->Suc) != Alg->FirstNode);
    for (i = *PopulationSize - 1; i >= 0 && Fitness.get()[i] > Cost; i--) {
		if ((d = DistanceToIndividual(i,Alg)) < MinDist) {
            CMin = i;
            MinDist = d;
        }
    }
    if (CMin == *PopulationSize - 1)
        return CMin;
    P = Population.get()[CMin];
    for (j = 0; j < Alg->Dimension; j++)
        Alg->NodeSet[P[j]].OldSuc = &Alg->NodeSet[P[j + 1]];
    for (i = 0; i < *PopulationSize; i++)
		if (i != CMin && (d = DistanceToIndividual(i,Alg)) < MinDist)
            return *PopulationSize - 1;
    return CMin;
}
