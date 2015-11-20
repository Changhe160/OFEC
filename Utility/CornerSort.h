//==================================================================
// This program was used in the following paper:
// H. Wang and X. Yao, ``Corner Sort for Pareto-based Many-Objective
// Optimization,'' IEEE Transactions on Cybernetics, 44(1):92-102, January
// 2014.
//
//==================================================================
#ifndef CONERSORT_H
#define CONERSORT_H

#include<malloc.h>
#include<stdlib.h>

namespace CornerSort{
struct node
{
	unsigned int index;//index in POP
	node *marknext;
	node *markpre;
	node *ranknext;
	node *rankpre;
};

//==================================================================
//unsigned int min_obj(double **POP,unsigned int obj_index,node **head,node **mh)
//POP-population??obj_index-objective index??head-linked list head?? mh-the head of the mark linked list
//Output:cur-the index of individual with the minimal obj_index-th objective
//==================================================================
unsigned int min_obj(double **POP,unsigned int obj_index,node **head,node **mh,int *comp)
{
	unsigned int i=0,cur;
	node *p=*mh,*q=NULL;
	cur=p->index;
	q=p;
	while(p->marknext!=NULL)
	{
		if (comp[obj_index] == 0) // minimization 
		{
			if (POP[p->marknext->index][obj_index] < POP[cur][obj_index])
			{
				cur = p->marknext->index;
				q = p->marknext;
			}
		}
		else{
			if (POP[p->marknext->index][obj_index] > POP[cur][obj_index])
			{
				cur = p->marknext->index;
				q = p->marknext;
			}
		}

		p=p->marknext;
	}
	if(q->markpre==NULL)
	{
		*mh=(*mh)->marknext;
		if((*mh)!=NULL)
		{
			(*mh)->markpre=NULL;
		}
	}
	else if(q->marknext==NULL)
	{
		q->markpre->marknext=NULL;
	}
	else
	{
		q->markpre->marknext=q->marknext;
		q->marknext->markpre=q->markpre;
	}
	if(q->rankpre==NULL)
	{
		*head=(*head)->ranknext;
		if((*head)!=NULL)
		{
			(*head)->rankpre=NULL;
		}
		free(q);
	}
	else if(q->ranknext==NULL)
	{
		p=q;
		q->rankpre->ranknext=NULL;
		free(p);
	}
	else
	{
		p=q;
		q->rankpre->ranknext=q->ranknext;
		q->ranknext->rankpre=q->rankpre;
		free(p);
	}
	return(cur);
}
//==================================================================
//void  cornerSort(double **POP,unsigned int m,unsigned int n,unsigned int *rank)
//POP-population??m-No. of Objectives??n-size of POP??rank-the ranking result as output
//==================================================================
void  cornerSort(double **POP, unsigned int m, unsigned int n, unsigned int *rank, int* comp)
{
	unsigned int i, j, cout = 0, markcout = 0, r = 0, cur, obj_index, flag;
	node *head = NULL, *mh = NULL, *p = NULL, *newone = NULL;
	//initialize linked list(rank and mark)
	for (i = 0; i < n; i++)
	{
		newone = (node*)malloc(sizeof(node));
		newone->index = i;
		newone->ranknext = NULL;
		if (i == 0)
		{
			head = newone;
			p = head;
			p->rankpre = NULL;
		}
		else
		{
			p->ranknext = newone;
			newone->rankpre = p;
			p = newone;
		}
	}

	while (head != NULL)
	{
		r = r + 1;//arrange the r-th ranking
		//unmark all the unranked solutions
		p = head;
		while (p != NULL)
		{
			p->marknext = p->ranknext;
			p->markpre = p->rankpre;
			p = p->ranknext;
		}
		mh = head;
		while (mh != NULL)//until all the solutions are marked
		{
			for (obj_index = 0; mh != NULL&&obj_index < m; ++obj_index)
			{
				// find solution of the best objective obj_index among unmarked ones
				cur = min_obj(POP, obj_index, &head, &mh, comp);// delete it in both mark and rank linked lists
				rank[cur] = r;
				cout++;
				p = mh;
				//mark the dominated solutions by cur
				while (p != NULL)
				{
					flag = 0;
					for (j = 0; j < m; ++j)
					{
						if (comp[j] == 0){
							if (j != obj_index&&POP[cur][j] > POP[p->index][j])
							{
								flag = 1;
								break;
							}
						}
						else{
							if (j != obj_index&&POP[cur][j] < POP[p->index][j])
							{
								flag = 1;
								break;
							}
						}

					}
					if (flag == 0)
					{
						if (p->markpre == NULL)
						{
							mh = mh->marknext;
							if (mh != NULL)
							{
								mh->markpre = NULL;
							}
						}
						else if (p->marknext == NULL)
						{
							p->markpre->marknext = NULL;
						}
						else
						{
							p->markpre->marknext = p->marknext;
							p->marknext->markpre = p->markpre;
						}
					}
					p = p->marknext;
				}
			}
		}
	}
}
}
#endif