#include "Nrand1DEIndividual.h"
#include "../../../Global/global.h"
Nrand1DEIndividual::Nrand1DEIndividual() :DEIndividual()
{
}

void Nrand1DEIndividual::recombine(double CR)
{
	int n,L,numDim;
	numDim=GET_NUM_DIM;
	n=Global::msp_global->getRandInt(0,numDim);
	L=0;
	do {
		L=L+1;
	}while(Global::msp_global->mp_uniformAlg->Next()<CR&&L<numDim);
	for(int i=0;i<numDim;i++)
	{
		if(i<L)
			m_pu.data()[(n+i)%numDim]=m_pv.data()[(n+i)%numDim];
		else
			m_pu.data()[(n+i)%numDim]=data()[(n+i)%numDim];
	}
}