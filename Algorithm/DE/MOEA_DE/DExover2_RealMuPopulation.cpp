#include "DExover2_RealMuPopulation.h"

void Diff_evo_xover2( const DEIndividual &parent1,const DEIndividual &parent2,const DEIndividual &parent3, DEIndividual &child);
void Real_mutation(DEIndividual &child,double r,double etam);

DExover2_RealMuPopulation::DExover2_RealMuPopulation(): DEPopulation<CodeVReal,DEIndividual>(),m_etam(20)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	m_r=1.0/numDim;
}

DExover2_RealMuPopulation::DExover2_RealMuPopulation(int popsize,bool mode): DEPopulation<CodeVReal,DEIndividual>(popsize,mode),m_etam(20)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	m_r=1.0/numDim;
}

void DExover2_RealMuPopulation::cross_mutate(const vector<int> &index, DEIndividual &child)
{
	Diff_evo_xover2(*m_pop[index[0]],*m_pop[index[1]],*m_pop[index[2]],child);
	Real_mutation(child,m_r,m_etam);
	child.evaluate();
}


void Diff_evo_xover2( const DEIndividual &parent1,const DEIndividual &parent2,const DEIndividual &parent3, DEIndividual &child)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	int idx_rnd = Global::msp_global->getRandInt(0,numDim);
	Boundary boun=CAST_PROBLEM_CONT->getSearchRange();
	double rate = 0.5;

	for(int n=0;n<numDim;n++)
	{
		/*Selected Two Parents*/
		child.data().m_x[n] = parent1.data().m_x[n] + rate*(parent3.data().m_x[n] - parent2.data().m_x[n]);

		if(child.data().m_x[n]<boun[n].m_lower){
			double rnd = Global::msp_global->mp_uniformAlg->Next();
			child.data().m_x[n] = boun[n].m_lower + rnd*(parent1.data().m_x[n] - boun[n].m_lower);
		}
		if(child.data().m_x[n]>boun[n].m_upper){ 
			double rnd = Global::msp_global->mp_uniformAlg->Next();
			child.data().m_x[n] = boun[n].m_upper - rnd*(boun[n].m_upper - parent1.data().m_x[n]);
		}	 
	}

}

void Real_mutation(DEIndividual &child,double r,double etam)
{
	int numDim=Global::msp_global->mp_problem->getNumDim();
	Boundary boun=CAST_PROBLEM_CONT->getSearchRange();

    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
	double eta_m = etam;

	int  id_rnd = Global::msp_global->getRandInt(0,numDim);

    for (int j=0; j<numDim; j++)
    {
		if (Global::msp_global->mp_uniformAlg->Next()<=r)
        {
            y  = child.data().m_x[j];
			yl = boun[j].m_lower;
			yu = boun[j].m_upper;
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd =Global::msp_global->mp_uniformAlg->Next();
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            child.data().m_x[j] = y;
        }
    }
}