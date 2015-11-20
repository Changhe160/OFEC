/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 21 September 2011
// Last modified: 12 Dec. 2014

#ifndef DEFINITION_H
#define DEFINITION_H
#include "include.h"

typedef map<string, int> STRING2ID;
typedef pair<string, string> ALG2PRO;
#define OFEC_PI acos(-1.0) ///3.14159265358979323846
#define OFEC_E exp(1.0)  // 2.71828182845904523536

#define MAX_NUM_RUN (Global::g_arg[param_numRun])
#define OFEC_CONSOLE
//#define OFEC_DEMON
//#define OFEC_DEBUG_ 
//#define OFEC_PROBLEM_DEBUG

enum Compare{MIN_OPT=0,MAX_OPT};
enum ProTag{SOP,MOP,DOP,MMP,SMP,ROOT,CONT,COMB,TSP,COP,VRP,TTP,JSP,KOP,SAT};
//SOP: single objective problem
//MOP: multi-objective problem
//DOP: dynamic optimization problem
//MMP: multi-modal problem
//SMP: single(uni) modal problem
//ROOT: robust optimzation problem
//CONT: continuous optimization problem
//COMB: combinatorial optimization problem
//TSP: travelling salesman problem
//COP: constraint optimization problem
//VRP: vehicle routing problem
//TTP: timetabling problem
//JSP: job shop problem
//KOP: knapsack optimization problem
//SAT: boolean satisfiability problem
enum ReturnFlag{Return_Normal=0,Return_Change,Return_Terminate,Return_ChangeNextEval,Return_Change_Timelinkage,Return_Loop\
	,Return_Change_Dim,Return_Error};
enum CompareResultFlag{Compare_Equal=0,Compare_Better,Compare_Worse,Compare_Non_Dominated,Compare_Non_Comparable};

enum PopInitMethod{POP_INIT_UNIFORM=0,POP_INIT_ORTHONORM,POP_INIT_CENTER,POP_INIT_USER_DEFINED,POP_INIT_RANDOM,POP_INIT_HEURIS};
enum ConvgProgeMode{PROGR_MEAN=0, PROGR_MEDEAN,PROGR_WORST, PROGR_BEST};
// the mode of the convergence speed graph,e.g., PROGR_MEAN denotes the mean value of all runs

enum SolutionValidation{VALIDATION_IGNORE=0,VALIDATION_REINITIALIZE,VALIDATION_REMAP,VALIDATION_SETTOBOUND};
//means of handling a solution that is out of the search space in continuous space

enum GAMutationStrategy{MUTAT_POLYNOMIAL,MUTAT_NORMAL,MUTAT_COMB};
enum GASelectionStrategy{SEL_TOURNAMENT=0,SEL_ROULETTE_WHEEL};
enum GAXoverStrategy{XOVER_ARITHMETIC,XOVER_SINGLEPOINT};

enum DEMutationStratgy{DE_rand_1, DE_best_1,DE_targetToBest_1,DE_best_2,DE_rand_2,DE_randToBest_1,DE_targetToRand_1};
enum DistanceMode{DIS_EUCLIDEAN=0,DIS_MANHATTAN,DIS_HAMMING};

enum MigrationMode{Migration_Best=0,Migration_Worst,Migration_Random};
enum MigrationTopology{Migration_Ring=0,Migration_Broadcast};
enum Param{param_numDim,param_numPeak,param_proName,param_algName,param_maxEvals,param_shiftLength,param_changeType,param_proId,\
	param_runId,param_algId,param_changeRatio,param_flagNumPeakChange,param_flagNumDimChange,param_peakNumChangeMode,param_flagNoise,\
	param_flagTimeLinkage,param_comDBGFunID,param_changeFre,param_noiseSeverity,param_timelinkageSeverity,param_popSize,param_evalCountFlag,\
	param_stepIndi,param_trainingTime,param_subPopSize,param_overlapDgre,param_convThreshold,param_timeWindow,param_sampleFre,param_gOptFlag,\
	param_workingDir,param_convFactor,param_numRun,param_numTask,param_exlRadius,param_solutionValidationMode,param_populationInitialMethod,\
	param_hibernatingRadius,param_minNumPopSize,param_numChange,param_numObj,param_case,param_resource4BestPop,param_proFileName,\
	param_xoverProbability, param_mutProbability, param_proTag, param_numGOpt, param_noiseFlag, param_predicFlag, param_changeType2, \
	param_peaksPerBox,param_interTest1, param_interTest2, param_interTest3, param_interTest4,param_numBox,param_heightConfigMode,param_peakCenter,\
	param_numParetoRegion,param_validRadius,param_attraction,param_radius,param_jumpHeight,param_variableRelation,param_peakShape,param_divisionMode,\
	param_peakOffset,param_flagIrregular,param_flagAsymmetric
};

enum ProgramMode{Program_Algorithm=0,Program_Problem};
//***********************Enviroment change types in GDBG system ******************************//
enum ChangeType{CT_SmallStep=0, CT_LargeStep,CT_Random,CT_Recurrent,CT_Chaotic,CT_RecurrentNoisy};
enum ComDBGFuncID{COMDBG_SPHERE=0, COMDBG_RASTRIGIN, COMDBG_GRIEWANK, COMDBG_ACKLEY, COMDBG_HYBRID};

typedef boost::unique_lock<boost::mutex> Ulock;
	
#define IS_PROBLEM_NAME(id,name)  (gGetProblemName(id).compare(name)==0) 
#define IS_ALG_NAME(id,name) (gGetAlgorithmName(id).compare(name)==0) 

#endif
