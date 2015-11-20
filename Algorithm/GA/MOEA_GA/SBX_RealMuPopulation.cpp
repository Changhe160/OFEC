#include "SBX_RealMuPopulation.h"
#include "../../../Utility/include.h"
#include "../../MOEA/NSGAIII/MathAux.h"



SBX_RealMuPopulation::SBX_RealMuPopulation(): GAPopulation<CodeVReal,GAIndividual<CodeVReal>>(),m_cr(1.0),m_ceta(30),m_meta(20)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	m_mr=1.0/numDim;
}

SBX_RealMuPopulation::SBX_RealMuPopulation(int popsize,bool mode): GAPopulation<CodeVReal,GAIndividual<CodeVReal>>(popsize,mode),\
	m_cr(1.0),m_ceta(30),m_meta(20)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	m_mr=1.0/numDim;
}


void SBX_RealMuPopulation::cross_mutate(const vector<int> &index, GAIndividual<CodeVReal> *child1, GAIndividual<CodeVReal> *child2)
{
	SimulatedBinaryCrossover(child1,child2,*m_pop[index[0]],*m_pop[index[1]],m_cr,m_ceta);

	PolynomialMutation(child1,m_mr,m_meta);
	PolynomialMutation(child2,m_mr,m_meta);

	child1->evaluate();
	child2->evaluate();
}


// ----------------------------------------------------------------------
// SimulatedBinaryCrossover : simulated binary crossover (SBX)
// The implementation was adapted from the code of function realcross() in crossover.c
// http://www.iitk.ac.in/kangal/codes/nsga2/nsga2-gnuplot-v1.1.6.tar.gz
//
// ref: http://www.slideshare.net/paskorn/simulated-binary-crossover-presentation#
// ----------------------------------------------------------------------
double get_betaq(double rand, double alpha, double ceta)
{
	double betaq = 0.0;
	if (rand <= (1.0/alpha))
	{
		betaq = std::pow((rand*alpha),(1.0/(ceta+1.0)));
	}
	else
	{
		betaq = std::pow((1.0/(2.0 - rand*alpha)),(1.0/(ceta+1.0)));
	}
	return betaq;
}

// ----------------------------------------------------------------------
void SimulatedBinaryCrossover(GAIndividual<CodeVReal> *child1, GAIndividual<CodeVReal> *child2, const GAIndividual<CodeVReal> &parent1, const GAIndividual<CodeVReal> &parent2,  double cr, double ceta)
{
	*child1 = parent1;
	*child2 = parent2;

	if (Global::msp_global->mp_uniformAlg->Next() > cr) return ; // not crossovered

	vector<double> &c1=child1->data().m_x, &c2=child2->data().m_x;
	const vector<double> &p1 = parent1.data().m_x, &p2 = parent2.data().m_x;						

	for (size_t i=0; i<c1.size(); i+=1)
	{
		if (Global::msp_global->mp_uniformAlg->Next() > 0.5) continue; // these two variables are not crossovered
		if (std::fabs(static_cast<double>(p1[i])-static_cast<double>(p2[i])) <= 1.0e-14) continue; // two values are the same
		
		double y1 = std::min(p1[i], p2[i]),
			   y2 = std::max(p1[i], p2[i]);

		double lb = CAST_PROBLEM_CONT->getSearchRange().getDomain(i).m_lower,
			   ub =CAST_PROBLEM_CONT->getSearchRange().getDomain(i).m_upper;

		double rand = Global::msp_global->mp_uniformAlg->Next();

		// child 1
		double beta = 1.0 + (2.0*(y1-lb)/(y2-y1)),
			   alpha = 2.0 - std::pow(beta, -(ceta+1.0));
		double betaq = get_betaq(rand, alpha, ceta);
		
		c1[i] = 0.5*((y1+y2)-betaq*(y2-y1));

		// child 2
		beta = 1.0 + (2.0*(ub-y2)/(y2-y1));
		alpha = 2.0 - std::pow(beta, -(ceta+1.0));
		betaq = get_betaq(rand, alpha, ceta);

		c2[i] = 0.5*((y1+y2)+betaq*(y2-y1));

		// boundary checking
		c1[i] = std::min(ub, std::max(lb, static_cast<double>(c1[i])));
		c2[i] = std::min(ub, std::max(lb, static_cast<double>(c2[i])));

		if (Global::msp_global->mp_uniformAlg->Next() <= 0.5)
		{
			std::swap(c1[i], c2[i]);
		}
	}

}// CSimulatedBinaryCrossover


// ----------------------------------------------------------------------
// The implementation was adapted from the code of function real_mutate_ind() in mutation.c in
// http://www.iitk.ac.in/kangal/codes/nsga2/nsga2-gnuplot-v1.1.6.tar.gz
//
// ref: http://www.slideshare.net/paskorn/simulated-binary-crossover-presentation#
// ---------------------------------------------------------------------
void PolynomialMutation(GAIndividual<CodeVReal> *indv, double mr, double meta)
{
    //int j;
    //double rnd, delta1, delta2, mut_pow, deltaq;
    //double y, yl, yu, val, xy;

	vector<double> &x=indv->data().m_x;

    for (size_t i=0; i<x.size(); i+=1)
    {
		if (Global::msp_global->mp_uniformAlg->Next() <= mr)
        {

            double y = x[i],
			       lb = CAST_PROBLEM_CONT->getSearchRange().getDomain(i).m_lower,
				   ub = CAST_PROBLEM_CONT->getSearchRange().getDomain(i).m_upper;

            double delta1 = (y-lb)/(ub-lb),
                   delta2 = (ub-y)/(ub-lb);
            
			double mut_pow = 1.0/(meta+1.0);

			double rnd = Global::msp_global->mp_uniformAlg->Next(), deltaq = 0.0;
            if (rnd <= 0.5)
            {
                double xy = 1.0-delta1;
                double val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(meta+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                double xy = 1.0-delta2;
                double val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(meta+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }

            y = y + deltaq*(ub-lb);
			y = std::min(ub, std::max(lb, y));

            x[i] = y;
        }
    }

}// CPolynomialMutation