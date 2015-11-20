#ifndef TYPES_ALG_PRO_H
#define TYPES_ALG_PRO_H
#include "../Utility/TypeList/Typelist.h"
#include "test.h"

typedef  FAMF<FAMFParticle,Swarm<CodeVReal,FAMFParticle>,FAMFPopPSO> ALG_FAMF_PSO;
typedef  FAMF<FAMFIndiDE,DEPopulation<CodeVReal,FAMFIndiDE>,FAMFPopDE> ALG_FAMF_DE;

typedef LOKI_TYPELIST_18(CPSOSwarm,CPSORSwarm,AMSOSwarm,SwarmLBest,NSGAII_DExover2RealMu,NSGAII_SBXRealMu,NSGAIII_DExover2RealMu,NSGAIII_SBXRealMu,MOEAD_DExover2RealMu
		,MOEAD_SBXRealMu,ALG_FAMF_PSO,ALG_FAMF_DE,SwarmGBest,SLPSO,AS,LKH::LKHAlg,DEBest2,DERand1) AlgList;
typedef LOKI_TYPELIST_50(FSphere,FNoncont_Rastrigin,FModified_Rastrigin,FRastrigin,FWeierstrass,FGriewank,FAckley,FStep,FQuartic_Noisy,FScaffer_F6,
		FRosenbrock,FSchwefel_2_13,FSchwefel_2_22,FSchwefel_1_2,FSchwefel_2_21,FSchwefel_2_6,FSchwefel,FPenalized_1,FPenalized_2,HybridComp,FElliptic,
		FMAX_global1,FMAX_global2,FMAX_global3,FMAX_global4,FMAX_global5,FGear_Train,FParEst_FMSoundWaves,FSix_humpCamelBack,FWaves,FBraninRCOS,
		FShubert,FMichalewicz,FValleys,FFive_hills,FCenter_peak,FKeane_Bump,FBeasley,FHimmenblau,FModified_Shekel,FSzu,FIBA,FShaffer,FVincent,
		CompositionDBG,RotationDBG,MovingPeak,ZDT1,ZDT2,ZDT3) ProList50;
typedef	Loki::TL::Append<ProList50,ZDT4>::Result ProList51;
typedef	Loki::TL::Append<ProList51,ZDT6>::Result ProList52;
typedef	Loki::TL::Append<ProList52,DTLZ1>::Result ProList53;
typedef	Loki::TL::Append<ProList53,DTLZ2>::Result ProList54;
typedef	Loki::TL::Append<ProList54,DTLZ3>::Result ProList55;
typedef	Loki::TL::Append<ProList55,DTLZ4>::Result ProList56;
typedef	Loki::TL::Append<ProList56,F1>::Result ProList57;
typedef	Loki::TL::Append<ProList57,F2>::Result ProList58;
typedef	Loki::TL::Append<ProList58,F3>::Result ProList59;
typedef	Loki::TL::Append<ProList59,F4>::Result ProList60;
typedef	Loki::TL::Append<ProList60,F5>::Result ProList61;
typedef	Loki::TL::Append<ProList61,F6>::Result ProList62;
typedef	Loki::TL::Append<ProList62,F7>::Result ProList63;
typedef	Loki::TL::Append<ProList63,F8>::Result ProList64;
typedef	Loki::TL::Append<ProList64,F9>::Result ProList65;
typedef	Loki::TL::Append<ProList65, TravellingSalesman>::Result ProList;
#endif

