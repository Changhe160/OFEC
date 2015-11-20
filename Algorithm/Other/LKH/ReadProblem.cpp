#include "LKH.h"
#include "Heap.h"

/*      
 * The ReadProblem function reads the problem data in TSPLIB format from the 
 * file specified in the parameter file (PROBLEM_FILE).
 *
 * The following description of the file format is extracted from the TSPLIB 
 * documentation.  
 *
 * The file consists of a specification part and a data part. The specification 
 * part contains information on the file format and on its contents. The data 
 * part contains explicit data.
 *
 * (1) The specification part
 *
 * All entries in this section are of the form <keyword> : <value>, where 
 * <keyword> denotes an alphanumerical keyword and <value> denotes 
 * alphanumerical or numerical data. The terms <string>, <integer> and <real> 
 * denote character string, integer or real data, respectively. The order of 
 * specification of the keywords in the data file is arbitrary (in principle), 
 * but must be consistent, i.e., whenever a keyword is specified, all 
 * necessary information for the correct interpretation of the keyword has to 
 * be known. 
 *
 * Below is given a list of all available keywords.
 *
 * NAME : <string>e
 * Identifies the data file.
 * 
 * TYPE : <string>
 * Specifies the type of data. Possible types are
 * TSP          Data for a symmetric traveling salesman problem
 * ATSP         Data for an asymmetric traveling salesman problem
 * HCP          Hamiltonian cycle problem data.
 * HPP          Hamiltonian path problem data (not available in TSPLIB)
 *
 * COMMENT : <string>
 * Additional comments (usually the name of the contributor or the creator of 
 * the problem instance is given here).
 *
 * DIMENSION : < integer>
 * The number of nodes.
 *
 * EDGE_WEIGHT_TYPE : <string>
 * Specifies how the edge weights (or distances) are given. The values are:
 * ATT          Special distance function for problem att48 and att532
 * CEIL_2D      Weights are Euclidean distances in 2-D rounded up
 * CEIL_3D      Weights are Euclidean distances in 3-D rounded up
 * EUC_2D       Weights are Euclidean distances in 2-D
 * EUC_3D       Weights are Euclidean distances in 3-D
 * EXPLICIT     Weights are listed explicitly in the corresponding section
 * GEO          Weights are geographical distances in kilometers (TSPLIB).
 *              Coordinates are given in the form DDD.MM where DDD are the
 *              degrees and MM the minutes
 * GEOM         Weights are geographical distances in meters (used for the 
 *              world TSP). Coordinates are given in decimal form    
 * GEO_MEEUS    Weights are geographical distances in kilometers, computed 
 *              according to Meeus' formula.  Coordinates are given in the 
 *              form DDD.MM where DDD are the degrees and MM the minutes
 * GEOM_MEEUS   Weights are geographical distances, computed according to 
 *              Meeus' formula. Coordinates are given in decimal form
 * MAN_2D       Weights are Manhattan distances in 2-D
 * MAN_3D       Weights are Manhattan distances in 3-D
 * MAX_2D       Weights are maximum distances in 2-D 
 * MAX_3D       Weights are maximum distances in 3-D
 * SPECIAL      There is a spcial distance function implemented in 
 *              the Distance_SPECIAL function.
 *
 * EDGE-WEIGHT_FORMAT : <string>
 * Describes the format of the edge weights if they are given explicitly. 
 * The values are
 * FUNCTION         Weights are given by a function (see above)
 * FULL_MATRIX      Weights are given by a full matrix
 * UPPER_ROW        Upper triangular matrix 
 *                      (row-wise without diagonal entries)
 * LOWER_ROW        Lower triangular matrix 
 *                      (row-wise without diagonal entries)     
 * UPPER_DIAG_ROW   Upper triangular matrix 
 *                      (row-wise including diagonal entries)
 * LOWER_DIAG_ROW   Lower triangular matrix 
 *                      (row-wise including diagonal entries)
 * UPPER_COL        Upper triangular matrix 
 *                      (column-wise without diagonal entries)
 * LOWER_COL        Lower triangular matrix 
 *                      (column-wise without diagonal entries)  
 * UPPER_DIAG_COL   Upper triangular matrix 
 *                      (column-wise including diagonal entries)
 * LOWER_DIAG_COL   Lower triangular matrix 
 *                      (column-wise including diagonal entries)
 *
 * EDGE_DATA_FORMAT : <string>
 * Describes the format in which the edges of a graph are given, if the 
 * graph is not complete. The values are
 * EDGE_LIST    The graph is given by an edge list
 * ADJ_LIST     The graph is given by an adjacency list
 *
 * NODE_COORD_TYPE : <string>
 * Specifies whether the coordinates are associated with each node 
 * (which, for example may be used for either graphical display or 
 * distance computations.
 * The values are
 * TWOD_COORDS      Nodes are specified by coordinates in 2-D
 * THREED_COORDS    Nodes are specified by coordinates in 3-D
 * NO_COORDS        The nodes do not have associated coordinates
 * The default value is NO_COORDS. In the current implementation, however, 
 * the value has no significance.
 *
 * DISPLAY_DATA_TYPE : <string>
 * Specifies how a graphical display of the nodes can be obtained. 
 * The values are
 * COORD_DISPLAY    Display is generated from the node coordinates
 * TWOD_DISPLAY     Explicit coordinates in 2-D are given
 * BO_DISPLAY       No graphical display is possible
 *
 * The default value is COORD_DISPLAY if node coordinates are specifies and 
 * NO_DISPLAY otherwise. In the current implementation, however, the value 
 * has no significance.
 *
 * EOF
 * Terminates input data. The entry is optional.
 *
 * (2) The data part
 *
 * Depending on the choice of specifications some additional data may be 
 * required. These data are given corresponding data sections following the 
 * specification part. Each data section begins with the corresponding 
 * keyword. The length of the sectionis either explicitly known form the 
 * format specification, or the section is terminated by an appropriate 
 * end-of-section identifier.
 *
 * NODE_COORD_SECTION :
 * Node coordinates are given in this section. Each line is of the form
 *  
 *      <integer> <real> <real>
 *       
 * if NODE_COORD_TYPE is TWOD_COORDS, or
 * 
 *      <integer> <real> <real> <real>
 *       
 * if NODE_COORD_TYPE is THREED_COORDS. The integers give the number of the 
 * respective nodes. The real numbers are the associated coordinates.
 *
 * EDGE_DATA_SECTION :
 * Edges of the graph are specified in either of the two formats allowed in 
 * the EDGE_DATA_FORAT entry. If a type is EDGE_LIST, then the edges are given 
 * as a sequence of lines of the form
 *  
 *      <integer> <integer>
 *       
 * each entry giving the terminal nodes of some edge. The list is terminated 
 * by a -1. If the type is ADJ_LIST, the section consists of adjacency lists 
 * for nodes.
 * The adjacency list of a node x is specified as
 * 
 *      <integer> <integer> ... <integer> -1
 *       
 * where the first integer gives the number of node x and the following 
 * integers (terminated by -1) the numbers of the nodes adjacent to x. 
 * The list of adjacency lists are terminated by an additional -1.
 *
 * FIXED_EDGES_SECTION :
 * In this section, edges are listed that are required to appear in each 
 * solution to the problem. The edges to be fixed are given in the form 
 * (per line)
 * 
 *      <integer> <integer>
 *       
 * meaning that the edge (arc) from the first node to the second node has 
 * to be contained in a solution. This section is terminated by a -1.
 *
 * DISPLAY_DATA_SECTION :
 * If DISPLAY_DATA_TYPE is TWOD_DISPLAY, the 2-dimensional coordinates from 
 * which a display can be generated are given in the form (per line)
 *  
 *      <integer> <real> <real>
 *       
 * The integers specify the respective nodes and the real numbers give the 
 * associated coordinates. The contents of this section, however, has no 
 * significance in the current implementation.
 *
 * TOUR_SECTION :
 * A tour is specified in this section. The tour is given by a list of 
 * integers giving the sequence in which the nodes are visited in the tour. 
 * The tour is terminated by a -1. Note: In contrast to the TSPLIB format, 
 * only one tour can be given in this section. The tour is used to limit 
 * the search (the last edge to be excluded in a non-gainful move must not 
 * belong to the tour). In addition, the Alpha field of its edges is set to 
 * -1.
 *
 * EDGE_WEIGHT_SECTION :
 * The edge weights are given in the format specifies by the EDGE_WEIGHT_FORMAT 
 * entry. At present, all explicit data are integral and is given in one of the
 * (self-explanatory) matrix formats, with explicitly known lengths.
 */

static const char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
static void CheckSpecificationPart(LKH::LKHAlg *Alg);
static char *Copy(char *S);
static void CreateNodes(LKH::LKHAlg *Alg);
static void Read_DIMENSION(LKH::LKHAlg *Alg);
static void Read_DISPLAY_DATA_SECTION(LKH::LKHAlg *Alg);
static void Read_DISPLAY_DATA_TYPE(LKH::LKHAlg *Alg);
static void Read_EDGE_DATA_FORMAT(LKH::LKHAlg *Alg);
static void Read_EDGE_DATA_SECTION(LKH::LKHAlg *Alg);
static void Read_EDGE_WEIGHT_FORMAT(LKH::LKHAlg *Alg);
static void Read_EDGE_WEIGHT_SECTION(LKH::LKHAlg *Alg);
static void Read_EDGE_WEIGHT_TYPE(LKH::LKHAlg *Alg);
static void Read_FIXED_EDGES_SECTION(LKH::LKHAlg *Alg);
static void Read_NAME(LKH::LKHAlg *Alg);
static void Read_NODE_COORD_SECTION(LKH::LKHAlg *Alg);
static void Read_NODE_COORD_TYPE(LKH::LKHAlg *Alg);
static void Read_TOUR_SECTION(FILE ** File,LKH::LKHAlg *Alg);
static void Read_TYPE(LKH::LKHAlg *Alg);
static int TwoDWeightType(LKH::LKHAlg *Alg);
static int ThreeDWeightType(LKH::LKHAlg *Alg);


void LKH::LKHAlg::ReadProblem()
{
    int i, K;
    char *Line, *Keyword;

	ostringstream os;
	os<<Global::g_arg[param_workingDir]<<"Problem/Combination/TSP/data/"<<ProblemFileName<<".tsp";
	if (!(ProblemFile = fopen(os.str().c_str(), "r")))
        eprintf("Cannot open PROBLEM_FILE: \"%s\"", os.str().c_str());
  //  if (TraceLevel >= 1)
   //     printff("Reading PROBLEM_FILE: \"%s\" ... ", os.str().c_str());
   // FreeStructures();
    FirstNode = 0;
    WeightType = WeightFormat = ProblemType = -1;
    CoordType = NO_COORDS;
    Name = "Unnamed";
    Type = EdgeWeightType = EdgeWeightFormat = 0;
    EdgeDataFormat = NodeCoordType = DisplayDataType = 0;
    Distance = 0;
    C = 0;
    c = 0;
    while ((Line = ReadLine(ProblemFile))) {
		if (!(Keyword = gStrtok_r(Line, Delimiters,&savePtr)))
            continue;
        for (i = 0; i < (int) strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "COMMENT"));
        else if (!strcmp(Keyword, "DEMAND_SECTION"))
            eprintf("Not implemented: %s", Keyword);
        else if (!strcmp(Keyword, "DEPOT_SECTION"))
            eprintf("Not implemented: %s", Keyword);
        else if (!strcmp(Keyword, "DIMENSION"))
			Read_DIMENSION(this);
        else if (!strcmp(Keyword, "DISPLAY_DATA_SECTION"))
            Read_DISPLAY_DATA_SECTION(this);
        else if (!strcmp(Keyword, "DISPLAY_DATA_TYPE"))
            Read_DISPLAY_DATA_TYPE(this);
        else if (!strcmp(Keyword, "EDGE_DATA_FORMAT"))
            Read_EDGE_DATA_FORMAT(this);
        else if (!strcmp(Keyword, "EDGE_DATA_SECTION"))
            Read_EDGE_DATA_SECTION(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_FORMAT"))
            Read_EDGE_WEIGHT_FORMAT(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_SECTION"))
            Read_EDGE_WEIGHT_SECTION(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_TYPE"))
            Read_EDGE_WEIGHT_TYPE(this);
        else if (!strcmp(Keyword, "EOF"))
            break;
        else if (!strcmp(Keyword, "FIXED_EDGES_SECTION"))
            Read_FIXED_EDGES_SECTION(this);
        else if (!strcmp(Keyword, "NAME"))
            Read_NAME(this);
        else if (!strcmp(Keyword, "NODE_COORD_SECTION"))
            Read_NODE_COORD_SECTION(this);
        else if (!strcmp(Keyword, "NODE_COORD_TYPE"))
            Read_NODE_COORD_TYPE(this);
        else if (!strcmp(Keyword, "TOUR_SECTION"))
            Read_TOUR_SECTION(&ProblemFile,this);
        else if (!strcmp(Keyword, "TYPE"))
            Read_TYPE(this);
        else
            eprintf("Unknown keyword: %s", Keyword);
    }
    Swaps = 0;

    /* Adjust parameters */
	Optimum=CAST_TSP->getGOpt().getGloObj()[0];
    if (Seed == 0)
        Seed = (unsigned) time(0);
    if (Precision == 0)
        Precision = 100;
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
  //  if (Runs == 0)
  //      Runs = 10;
	if(Runs!=1)
		Runs=1;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (SubproblemSize >= Dimension)
        SubproblemSize = Dimension;
    else if (SubproblemSize == 0) {
        if (AscentCandidates > Dimension - 1)
            AscentCandidates = Dimension - 1;
        if (InitialPeriod < 0) {
            InitialPeriod = Dimension / 2;
            if (InitialPeriod < 100)
                InitialPeriod = 100;
        }
        if (Excess < 0)
            Excess = 1.0 / Dimension;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        MakeHeap(Dimension,this);
    }

    if (CostMatrix == 0 && Dimension <= MaxMatrixDimension && Distance != 0
        && Distance != &LKH::LKHAlg::Distance_1 && Distance != &LKH::LKHAlg::Distance_ATSP &&
        WeightType != GEO && WeightType != GEOM &&
        WeightType != GEO_MEEUS && WeightType != GEOM_MEEUS) {
        Node *Ni, *Nj;
        assert(CostMatrix =
               (int *) calloc((size_t) Dimension * (Dimension - 1) / 2,
                              sizeof(int)));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
            if (ProblemType != HPP || Ni->Id < Dimension)
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : (this->*Distance)(Ni, Nj);
            else
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = 0;
        }
        while ((Ni = Ni->Suc) != FirstNode);
        WeightType = EXPLICIT;
        c = 0;
    }
    if (Precision > 1 && (WeightType == EXPLICIT || ProblemType == ATSP)) {
        int j, n = ProblemType == ATSP ? Dimension / 2 : Dimension;
        for (i = 2; i <= n; i++) {
            Node *N = &NodeSet[i];
            for (j = 1; j < i; j++)
                if (N->C[j] * Precision / Precision != N->C[j])
                    eprintf("PRECISION (= %d) is too large", Precision);
        }
    }
    C = WeightType == EXPLICIT ? &LKH::LKHAlg::C_EXPLICIT : &LKH::LKHAlg::C_FUNCTION;
    D = WeightType == EXPLICIT ? &LKH::LKHAlg::D_EXPLICIT : &LKH::LKHAlg::D_FUNCTION;
    if (SubsequentMoveType == 0)
        SubsequentMoveType = MoveType;
    K = MoveType >= SubsequentMoveType
        || !SubsequentPatching ? MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1 && NonsequentialMoveType >= 4) {
        BestMove = BestSubsequentMove = &LKH::LKHAlg::BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
                { 0, 0, &LKH::LKHAlg::Best2OptMove, &LKH::LKHAlg::Best3OptMove,
                &LKH::LKHAlg::Best4OptMove, &LKH::LKHAlg::Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, &LKH::LKHAlg::Best2OptMove, &LKH::LKHAlg::Best3OptMove,
            &LKH::LKHAlg::Best4OptMove, &LKH::LKHAlg::Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : &LKH::LKHAlg::BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : &LKH::LKHAlg::BestKOptMove;
    }
    if (ProblemType == HCP || ProblemType == HPP)
        MaxCandidates = 0;
    if (TraceLevel >= 1) {
      //  printff("done\n");
       // PrintParameters();
    } else
        printff("PROBLEM_FILE = %s\n",
                ProblemFileName ? ProblemFileName : "");
    fclose(ProblemFile);
    if (InitialTourFileName)
        ReadTour(InitialTourFileName, &InitialTourFile);
    if (InputTourFileName)
        ReadTour(InputTourFileName, &InputTourFile);
    if (SubproblemTourFileName && SubproblemSize > 0)
        ReadTour(SubproblemTourFileName, &SubproblemTourFile);
    if (MergeTourFiles >= 1) {
        free(MergeTourFile);
        assert(MergeTourFile =
               (FILE **) malloc(MergeTourFiles * sizeof(FILE *)));
        for (i = 0; i < MergeTourFiles; i++)
            ReadTour(MergeTourFileName[i], &MergeTourFile[i]);
    }
    free(LastLine);
    LastLine = 0;
}

static int TwoDWeightType(LKH::LKHAlg *Alg)
{
    return Alg->WeightType == LKH::EUC_2D || Alg->WeightType == LKH::MAX_2D ||
        Alg->WeightType == LKH::MAN_2D || Alg->WeightType == LKH::CEIL_2D ||
        Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEOM ||
        Alg->WeightType == LKH::GEO_MEEUS || Alg->WeightType == LKH::GEOM_MEEUS ||
        Alg->WeightType == LKH::ATT ||
        (Alg->WeightType == LKH::SPECIAL && Alg->CoordType == LKH::TWOD_COORDS);
}

static int ThreeDWeightType(LKH::LKHAlg *Alg)
{
    return Alg->WeightType == LKH::EUC_3D || Alg->WeightType == LKH::MAX_3D ||
        Alg->WeightType == LKH::MAN_3D || Alg->WeightType == LKH::CEIL_3D ||
        (Alg->WeightType == LKH::SPECIAL && Alg->CoordType == LKH::THREED_COORDS);
}

static void CheckSpecificationPart(LKH::LKHAlg *Alg)
{
    if (Alg->ProblemType == -1)
        Alg->eprintf("TYPE is missing");
    if (Alg->Dimension < 3)
        Alg->eprintf("DIMENSION < 3 or not specified");
    if (Alg->WeightType == -1 && Alg->ProblemType != LKH::ATSP && Alg->ProblemType != LKH::HCP)
        Alg->eprintf("EDGE_WEIGHT_TYPE is missing");
    if (Alg->WeightType == LKH::EXPLICIT && Alg->WeightFormat == -1)
        Alg->eprintf("EDGE_WEIGHT_FORMAT is missing");
    if (Alg->WeightType == LKH::EXPLICIT && Alg->WeightFormat == LKH::FUNCTION)
        Alg->eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (Alg->WeightType != LKH::EXPLICIT
        && (Alg->WeightType != LKH::SPECIAL || Alg->CoordType != LKH::NO_COORDS)
        && Alg->WeightType != -1 && Alg->WeightFormat != -1
        && Alg->WeightFormat != LKH::FUNCTION)
        Alg->eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (Alg->ProblemType == LKH::ATSP && Alg->WeightType != LKH::EXPLICIT && Alg->WeightType != -1)
        Alg->eprintf("Conflicting TYPE and EDGE_WEIGHT_TYPE");
    if (Alg->ProblemType == LKH::ATSP && Alg->WeightFormat != LKH::FULL_MATRIX)
        Alg->eprintf("Conflicting TYPE and EDGE_WEIGHT_FORMAT");
    if (Alg->CandidateSetType == LKH::DELAUNAY && !TwoDWeightType(Alg)
        && Alg->MaxCandidates > 0)
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = DELAUNAY");
    if (Alg->CandidateSetType == LKH::NN && !TwoDWeightType(Alg)
        && !ThreeDWeightType(Alg) && Alg->MaxCandidates > 0)
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = "
             "NEAREST-NEIGHBOR");
    if (Alg->CandidateSetType == LKH::QUADRANT && !TwoDWeightType(Alg)
        && !ThreeDWeightType(Alg) && Alg->MaxCandidates + Alg->ExtraCandidates > 0)
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = QUADRANT");
    if (Alg->ExtraCandidateSetType == LKH::NN && !TwoDWeightType(Alg)
        && !ThreeDWeightType(Alg) && Alg->ExtraCandidates > 0)
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "NEAREST-NEIGHBOR");
    if (Alg->ExtraCandidateSetType == LKH::QUADRANT && !TwoDWeightType(Alg)
        && !ThreeDWeightType(Alg)
        && Alg->ExtraCandidates > 0)
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "QUADRANT");
	if (Alg->InitialTourAlgorithm == LKH::QUICK_BORUVKA && !TwoDWeightType(Alg)
		&& !ThreeDWeightType(Alg))
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "QUICK-BORUVKA");
	if (Alg->InitialTourAlgorithm == LKH::SIERPINSKI && !TwoDWeightType(Alg))
        Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "SIERPINSKI");
	if (Alg->DelaunayPartitioning && !TwoDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for DELAUNAY specification");
	if (Alg->KarpPartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for KARP specification");
	if (Alg->KMeansPartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for K-MEANS specification");
	if (Alg->MoorePartitioning && !TwoDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for MOORE specification");
	if (Alg->RohePartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for ROHE specification");
	if (Alg->SierpinskiPartitioning && !TwoDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for SIERPINSKI specification");
	if (Alg->SubproblemBorders && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
        Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for BORDERS specification");
}

static char *Copy(char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    assert(Buffer = (char *) malloc(strlen(S) + 1));
    strcpy(Buffer, S);
    return Buffer;
}

static void CreateNodes(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Prev = 0, *N = 0;
    int i;

    if (Alg->Dimension <= 0)
        Alg->eprintf("DIMENSION is not positive (or not specified)");
    if (Alg->ProblemType == LKH::ATSP)
        Alg->Dimension *= 2;
    else if (Alg->ProblemType == LKH::HPP) {
        Alg->Dimension++;
        if (Alg->Dimension >Alg->MaxMatrixDimension)
            Alg->eprintf("Dimension too large in HPP problem");
    }
    assert(Alg->NodeSet = (LKH::LKHAlg::Node *) calloc(Alg->Dimension + 1, sizeof(LKH::LKHAlg::Node)));
    for (i = 1; i <= Alg->Dimension; i++, Prev = N) {
        N = &Alg->NodeSet[i];
        if (i == 1)
            Alg->FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
        if (Alg->MergeTourFiles >= 1)
            assert(N->MergeSuc =
                   (LKH::LKHAlg::Node **) calloc(Alg->MergeTourFiles, sizeof(LKH::LKHAlg::Node *)));
    }
    Link(N, Alg->FirstNode);
}

static void Read_NAME(LKH::LKHAlg *Alg)
{
    if (!(Alg->Name = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("NAME: string expected");
}

static void Read_DIMENSION(LKH::LKHAlg *Alg)
{
    char *Token = gStrtok_r(0, Delimiters,&Alg->savePtr);

    if (!Token || !sscanf(Token, "%d", &Alg->Dimension))
        Alg->eprintf("DIMENSION: integer expected");
    Alg->DimensionSaved = Alg->Dimension;
}

static void Read_DISPLAY_DATA_SECTION(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *N;
    int Id, i;

	CheckSpecificationPart(Alg);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (!Alg->DisplayDataType || strcmp(Alg->DisplayDataType, "TWOD_DISPLAY"))
        Alg->eprintf
            ("DISPLAY_DATA_SECTION conflicts with DISPLAY_DATA_TYPE: %s",
             Alg->DisplayDataType);
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    N = Alg->FirstNode;
    for (i = 1; i <= Alg->Dimension; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
            Alg->eprintf("Missing nodes in DIPLAY_DATA_SECTION");
        if (Id <= 0 || Id >Alg-> Dimension)
           Alg-> eprintf("(DIPLAY_DATA_SECTION) Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
            Alg->eprintf("(DIPLAY_DATA_SECTION) Node number occours twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%lf", &N->X))
            Alg->eprintf("Missing X-coordinate in DIPLAY_DATA_SECTION");
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Y))
            Alg->eprintf("Missing Y-coordinate in DIPLAY_DATA_SECTION");
    }
    N = Alg->FirstNode;
    do
        if (!N->V)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
        Alg->eprintf("(DIPLAY_DATA_SECTION) No coordinates given for node %d",
                N->Id);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
}

static void Read_DISPLAY_DATA_TYPE(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->DisplayDataType = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("DISPLAY_DATA_TYPE: string expected");
    for (i = 0; i < strlen(Alg->DisplayDataType); i++)
        Alg->DisplayDataType[i] = (char) toupper(Alg->DisplayDataType[i]);
    if (strcmp(Alg->DisplayDataType, "COORD_DISPLAY") &&
        strcmp(Alg->DisplayDataType, "TWOD_DISPLAY") &&
        strcmp(Alg->DisplayDataType, "NO_DISPLAY"))
        Alg->eprintf("Unknown DISPLAY_DATA_TYPE: %s", Alg->DisplayDataType);
}

static void Read_EDGE_DATA_FORMAT(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeDataFormat = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("EDGE_DATA_FORMAT: string expected");
    for (i = 0; i < strlen(Alg->EdgeDataFormat); i++)
        Alg->EdgeDataFormat[i] = (char) toupper(Alg->EdgeDataFormat[i]);
    if (strcmp(Alg->EdgeDataFormat, "EDGE_LIST") &&
        strcmp(Alg->EdgeDataFormat, "ADJ_LIST"))
        Alg->eprintf("Unknown EDGE_DATA_FORMAT: %s", Alg->EdgeDataFormat);
    if (Alg->SubproblemTourFileName)
        Alg->eprintf("(EDGE_DATA_FORMAT)"
                " cannot be used together with SUBPROBLEM_TOUR_FILE");
}

static void Read_EDGE_DATA_SECTION(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Ni, *Nj;
    int i, j;

	CheckSpecificationPart(Alg);
    if (!Alg->EdgeDataFormat)
        Alg->eprintf("Missing EDGE_DATA_FORMAT specification");
    if (!Alg->FirstNode)
		CreateNodes(Alg);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (!strcmp(Alg->EdgeDataFormat, "EDGE_LIST")) {
        if (!Alg->fscanint(Alg->ProblemFile, &i))
            i = -1;
        while (i != -1) {
            if (i <= 0
                || i > (Alg->ProblemType != LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
                Alg->eprintf("(EDGE_DATA_SECTION) Node number out of range: %d",
                        i);
            Alg->fscanint(Alg->ProblemFile, &j);
            if (j <= 0
                || j > ( Alg->ProblemType != LKH::ATSP ?  Alg->Dimension :  Alg->Dimension / 2))
                 Alg->eprintf("(EDGE_DATA_SECTION) Node number out of range: %d",
                        j);
            if (i == j)
                 Alg->eprintf("(EDGE_DATA_SECTION) Illgal edge: %d to %d", i, j);
            if ( Alg->ProblemType == LKH::ATSP) {
                i +=  Alg->Dimension / 2;
                j +=  Alg->Dimension / 2;
            }
            Ni = & Alg->NodeSet[i];
            Nj = & Alg->NodeSet[j];
            if (!Ni->CandidateSet) {
                assert(Ni->CandidateSet =
                       (LKH::LKHAlg::Candidate *) calloc(3, sizeof(LKH::LKHAlg::Candidate)));
                Ni->CandidateSet[0].To = Nj;
                Ni->CandidateSet[0].Cost = 0;
                Ni->CandidateSet[0].Alpha = 0;
                Ni->V = 1;
            } else if (!Alg->IsCandidate(Ni, Nj)) {
                Ni->CandidateSet[Ni->V].To = Nj;
                Ni->CandidateSet[Ni->V].Cost = 0;
                Ni->CandidateSet[Ni->V].Alpha = Ni->V;
                assert(Ni->CandidateSet =
                       (LKH::LKHAlg::Candidate *) realloc(Ni->CandidateSet,
                                             (++Ni->V +
                                              1) * sizeof(LKH::LKHAlg::Candidate)));
                Ni->CandidateSet[Ni->V].To = 0;
            }
            if (Alg->ProblemType != LKH::ATSP) {
                if (!Nj->CandidateSet) {
                    assert(Nj->CandidateSet =
                           (LKH::LKHAlg::Candidate *) calloc(3, sizeof(LKH::LKHAlg::Candidate)));
                    Nj->CandidateSet[0].To = Ni;
                    Nj->CandidateSet[0].Cost = 0;
                    Nj->CandidateSet[0].Alpha = 0;
                    Nj->V = 1;
                } else if (Alg->ProblemType !=!Alg->IsCandidate(Nj, Ni)) {
                    Nj->CandidateSet[Nj->V].To = Ni;
                    Nj->CandidateSet[Nj->V].Cost = 0;
                    Nj->CandidateSet[Nj->V].Alpha = Nj->V;
                    assert(Nj->CandidateSet =
                           (LKH::LKHAlg::Candidate *) realloc(Nj->CandidateSet,
                                                 (++Nj->V +
                                                  1) * sizeof(LKH::LKHAlg::Candidate)));
                    Nj->CandidateSet[Nj->V].To = 0;
                }
            }
            if (!Alg->fscanint(Alg->ProblemFile, &i))
                i = -1;
        }
    } else if (!strcmp(Alg->EdgeDataFormat, "ADJ_LIST")) {
        Ni = Alg->FirstNode;
        do
            Ni->V = 0;
        while ((Ni = Ni->Suc) != Alg->FirstNode);
        if (!Alg->fscanint(Alg->ProblemFile, &i))
            i = -1;
        while (i != -1) {
            if (i <= 0
                || i > (Alg->ProblemType != LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
                Alg->eprintf("(EDGE_DATA_SECTION) Node number out of range: %d",
                        i);
            if (Alg->ProblemType == LKH::ATSP)
                i += Alg->Dimension / 2;
            Ni = &Alg->NodeSet[i];
            Alg->fscanint(Alg->ProblemFile, &j);
            while (j != -1) {
                if (j <= 0
                    || j > (Alg->ProblemType !=
                            LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
                   Alg-> eprintf
                        ("(EDGE_DATA_SECTION) Node number out of range: %d",
                         j);
                if (i == j)
                    Alg->eprintf("(EDGE_DATA_SECTION) Illgal edge: %d to %d",
                            i, j);
                if (Alg->ProblemType == LKH::ATSP)
                    j += Alg->Dimension / 2;
                Nj = &Alg->NodeSet[j];
                if (!Ni->CandidateSet) {
                    assert(Ni->CandidateSet =
                           (LKH::LKHAlg::Candidate *) calloc(3, sizeof(LKH::LKHAlg::Candidate)));
                    Ni->CandidateSet[0].To = Nj;
                    Ni->CandidateSet[0].Cost = 0;
                    Ni->CandidateSet[0].Alpha = 0;
                    Ni->V = 1;
                } else if (!Alg->IsCandidate(Ni, Nj)) {
                    Ni->CandidateSet[Ni->V].To = Nj;
                    Ni->CandidateSet[Ni->V].Cost = 0;
                    Ni->CandidateSet[Ni->V].Alpha = Ni->V;
                    assert(Ni->CandidateSet =
                           (LKH::LKHAlg::Candidate *) realloc(Ni->CandidateSet,
                                                 (++Ni->V +
                                                  1) * sizeof(LKH::LKHAlg::Candidate)));
                    Ni->CandidateSet[Ni->V].To = 0;
                }
                if (Alg->ProblemType != LKH::ATSP) {
                    if (!Nj->CandidateSet) {
                        assert(Nj->CandidateSet =
                               (LKH::LKHAlg::Candidate *) calloc(3, sizeof(LKH::LKHAlg::Candidate)));
                        Nj->CandidateSet[0].To = Ni;
                        Nj->CandidateSet[0].Cost = 0;
                        Nj->CandidateSet[0].Alpha = 0;
                        Nj->V = 1;
                    } else if (!Alg->IsCandidate(Nj, Ni)) {
                        Nj->CandidateSet[Nj->V].To = Ni;
                        Nj->CandidateSet[Nj->V].Cost = 0;
                        Nj->CandidateSet[Nj->V].Alpha = Nj->V;
                        assert(Nj->CandidateSet =
                               (LKH::LKHAlg::Candidate *) realloc(Nj->CandidateSet,
                                                     (++Nj->V +
                                                      1) * sizeof(LKH::LKHAlg::Candidate)));
                        Nj->CandidateSet[Nj->V].To = 0;
                    }
                }
                if (!Alg->fscanint(Alg->ProblemFile, &j))
                    j = -1;
            }
            if (!Alg->fscanint(Alg->ProblemFile, &i))
                i = -1;
        }
    } else
        Alg->eprintf("(EDGE_DATA_SECTION) No EDGE_DATA_FORMAT specified");
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
	Alg->Distance = &LKH::LKHAlg::Distance_1;
}

static void Read_EDGE_WEIGHT_FORMAT(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeWeightFormat = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("EDGE_WEIGHT_FORMAT: string expected");
    for (i = 0; i < strlen(Alg->EdgeWeightFormat); i++)
        Alg->EdgeWeightFormat[i] = (char) toupper(Alg->EdgeWeightFormat[i]);
    if (!strcmp(Alg->EdgeWeightFormat, "FUNCTION"))
        Alg->WeightFormat = LKH::FUNCTION;
    else if (!strcmp(Alg->EdgeWeightFormat, "FULL_MATRIX"))
        Alg->WeightFormat = LKH::FULL_MATRIX;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_ROW"))
       Alg-> WeightFormat = LKH::UPPER_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_ROW"))
        Alg->WeightFormat = LKH::LOWER_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_DIAG_ROW"))
        Alg->WeightFormat = LKH::UPPER_DIAG_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_DIAG_ROW"))
        Alg->WeightFormat = LKH::LOWER_DIAG_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_COL"))
        Alg->WeightFormat = LKH::UPPER_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_COL"))
        Alg->WeightFormat = LKH::LOWER_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_DIAG_COL"))
        Alg->WeightFormat = LKH::UPPER_DIAG_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_DIAG_COL"))
        Alg->WeightFormat = LKH::LOWER_DIAG_COL;
    else
        Alg->eprintf("Unknown EDGE_WEIGHT_FORMAT: %s", Alg->EdgeWeightFormat);
}

static void Read_EDGE_WEIGHT_SECTION(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Ni, *Nj;
    int i, j, n, W;

	CheckSpecificationPart(Alg);
    if (!Alg->FirstNode)
		CreateNodes(Alg);
    if (Alg->ProblemType != LKH::ATSP) {
        assert(Alg->CostMatrix =
               (int *) calloc((size_t) Alg->Dimension * (Alg->Dimension - 1) / 2,
                              sizeof(int)));
        Ni = Alg->FirstNode->Suc;
        do {
            Ni->C =
                &Alg->CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
        }
        while ((Ni = Ni->Suc) != Alg->FirstNode);
    } else {
        n = Alg->Dimension / 2;
        assert(Alg->CostMatrix = (int *) calloc((size_t) n * n, sizeof(int)));
        for (Ni = Alg->FirstNode; Ni->Id <= n; Ni = Ni->Suc)
            Ni->C = &Alg->CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    switch (Alg->WeightFormat) {
    case LKH::FULL_MATRIX:
        if (Alg->ProblemType == LKH::ATSP) {
            n = Alg->Dimension / 2;
            for (i = 1; i <= n; i++) {
                Ni = &Alg->NodeSet[i];
                for (j = 1; j <= n; j++) {
                    if (!Alg->fscanint(Alg->ProblemFile, &W))
                        Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    Ni->C[j] = W;
                    if (i != j && W > Alg->M)
                        Alg->M = W;
                }
                Nj = &Alg->NodeSet[i + n];
                if (!Ni->FixedTo1)
                    Ni->FixedTo1 = Nj;
                else if (!Ni->FixedTo2)
                    Ni->FixedTo2 = Nj;
                if (!Nj->FixedTo1)
                    Nj->FixedTo1 = Ni;
                else if (!Nj->FixedTo2)
                    Nj->FixedTo2 = Ni;
            }
            Alg->Distance = &LKH::LKHAlg::Distance_ATSP;
            Alg->WeightType = -1;
        } else
            for (i = 1, Ni = Alg->FirstNode; i <= Alg->Dimension; i++, Ni = Ni->Suc) {
                for (j = 1; j <= Alg->Dimension; j++) {
                    if (!Alg->fscanint(Alg->ProblemFile, &W))
                        Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                    if (j < i)
                        Ni->C[j] = W;
                }
            }
        break;
    case LKH::UPPER_ROW:
        for (i = 1, Ni = Alg->FirstNode; i < Alg->Dimension; i++, Ni = Ni->Suc) {
            for (j = i + 1, Nj = Ni->Suc; j <= Alg->Dimension;
                 j++, Nj = Nj->Suc) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                Nj->C[i] = W;
            }
        }
        break;
    case LKH::LOWER_ROW:
        for (i = 2, Ni = Alg->FirstNode->Suc; i <= Alg->Dimension; i++, Ni = Ni->Suc) {
            for (j = 1; j < i; j++) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                Ni->C[j] = W;
            }
        }
        break;
    case LKH::UPPER_DIAG_ROW:
        for (i = 1, Ni = Alg->FirstNode; i <= Alg->Dimension; i++, Ni = Ni->Suc) {
            for (j = i, Nj = Ni; j <= Alg->Dimension; j++, Nj = Nj->Suc) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                if (i != j)
                    Nj->C[i] = W;
            }
        }
        break;
    case LKH::LOWER_DIAG_ROW:
        for (i = 1, Ni = Alg->FirstNode; i <= Alg->Dimension; i++, Ni = Ni->Suc) {
            for (j = 1; j <= i; j++) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                if (j != i)
                    Ni->C[j] = W;
            }
        }
        break;
    case LKH::UPPER_COL:
        for (j = 2, Nj = Alg->FirstNode->Suc; j <= Alg->Dimension; j++, Nj = Nj->Suc) {
            for (i = 1; i < j; i++) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                Nj->C[i] = W;
            }
        }
        break;
    case LKH::LOWER_COL:
        for (j = 1, Nj = Alg->FirstNode; j < Alg->Dimension; j++, Nj = Nj->Suc) {
            for (i = j + 1, Ni = Nj->Suc; i <= Alg->Dimension;
                 i++, Ni = Ni->Suc) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                Ni->C[j] = W;
            }
        }
        break;
    case LKH::UPPER_DIAG_COL:
        for (j = 1, Nj = Alg->FirstNode; j <= Alg->Dimension; j++, Nj = Nj->Suc) {
            for (i = 1; i <= j; i++) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                if (i != j)
                    Nj->C[i] = W;
            }
        }
        break;
    case LKH::LOWER_DIAG_COL:
        for (j = 1, Nj = Alg->FirstNode; j <= Alg->Dimension; j++, Nj = Nj->Suc) {
            for (i = j, Ni = Nj; i <= Alg->Dimension; i++, Ni = Ni->Suc) {
                if (!Alg->fscanint(Alg->ProblemFile, &W))
                    Alg->eprintf("Missing weight in EDGE_WEIGHT_SECTION");
                if (i != j)
                    Ni->C[j] = W;
            }
        }
        break;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
}

static void Read_EDGE_WEIGHT_TYPE(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeWeightType = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("EDGE_WEIGHT_TYPE: string expected");
    for (i = 0; i < strlen(Alg->EdgeWeightType); i++)
        Alg->EdgeWeightType[i] = (char) toupper(Alg->EdgeWeightType[i]);
    if (!strcmp(Alg->EdgeWeightType, "ATT")) {
        Alg->WeightType = LKH::ATT;
        Alg->Distance = &LKH::LKHAlg::Distance_ATT;
        Alg->c = &LKH::LKHAlg::c_ATT;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "CEIL_2D")) {
        Alg->WeightType = LKH::CEIL_2D;
        Alg->Distance = &LKH::LKHAlg::Distance_CEIL_2D;
        Alg->c = &LKH::LKHAlg::c_CEIL_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "CEIL_3D")) {
        Alg->WeightType = LKH::CEIL_3D;
        Alg->Distance = &LKH::LKHAlg::Distance_CEIL_3D;
        Alg->c = &LKH::LKHAlg::c_CEIL_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "EUC_2D")) {
        Alg->WeightType = LKH::EUC_2D;
        Alg->Distance = &LKH::LKHAlg::Distance_EUC_2D;
        Alg->c = &LKH::LKHAlg::c_EUC_2D;
       Alg-> CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "EUC_3D")) {
       Alg-> WeightType = LKH::EUC_3D;
        Alg->Distance = &LKH::LKHAlg::Distance_EUC_3D;
        Alg->c = &LKH::LKHAlg::c_EUC_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "EXPLICIT")) {
        Alg->WeightType = LKH::EXPLICIT;
        Alg->Distance = &LKH::LKHAlg::Distance_EXPLICIT;
    } else if (!strcmp(Alg->EdgeWeightType, "MAN_2D")) {
        Alg->WeightType = LKH::MAN_2D;
        Alg->Distance = &LKH::LKHAlg::Distance_MAN_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAN_3D")) {
        Alg->WeightType = LKH::MAN_3D;
        Alg->Distance = &LKH::LKHAlg::Distance_MAN_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAX_2D")) {
       Alg-> WeightType = LKH::MAX_2D;
        Alg->Distance = &LKH::LKHAlg::Distance_MAX_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAX_3D")) {
        Alg->WeightType = LKH::MAX_3D;
        Alg->Distance = &LKH::LKHAlg::Distance_MAX_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEO")) {
        Alg->WeightType = LKH::GEO;
        Alg->Distance = &LKH::LKHAlg::Distance_GEO;
        Alg->c = &LKH::LKHAlg::c_GEO;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEOM")) {
        Alg->WeightType = LKH::GEOM;
        Alg->Distance = &LKH::LKHAlg::Distance_GEOM;
        Alg->c = &LKH::LKHAlg::c_GEOM;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEO_MEEUS")) {
        Alg->WeightType = LKH::GEO_MEEUS;
        Alg->Distance = &LKH::LKHAlg::Distance_GEO_MEEUS;
        Alg->c = &LKH::LKHAlg::c_GEO_MEEUS;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEOM_MEEUS")) {
        Alg->WeightType = LKH::GEOM_MEEUS;
        Alg->Distance = &LKH::LKHAlg::Distance_GEOM_MEEUS;
        Alg->c = &LKH::LKHAlg::c_GEOM_MEEUS;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "SPECIAL")) {
        Alg->WeightType = LKH::SPECIAL;
        Alg->Distance = &LKH::LKHAlg::Distance_SPECIAL;
    } else if (!strcmp(Alg->EdgeWeightType, "XRAY1") ||
               !strcmp(Alg->EdgeWeightType, "XRAY2"))
        Alg->eprintf("EDGE_WEIGHT_TYPE not implemented: %s", Alg->EdgeWeightType);
    else
        Alg->eprintf("Unknown EDGE_WEIGHT_TYPE: %s", Alg->EdgeWeightType);
}

static void Read_FIXED_EDGES_SECTION(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *Ni, *Nj, *N, *NPrev = 0, *NNext;
    int i, j, Count = 0;

	CheckSpecificationPart(Alg);
    if (!Alg->FirstNode)
		CreateNodes(Alg);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (!Alg->fscanint(Alg->ProblemFile, &i))
        i = -1;
    while (i != -1) {
        if (i <= 0
            || i > (Alg->ProblemType != LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
            Alg->eprintf("(FIXED_EDGES_SECTION) Node number out of range: %d",
                    i);
        Alg->fscanint(Alg->ProblemFile, &j);
        if (j <= 0
            || j > (Alg->ProblemType != LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
            Alg->eprintf("(FIXED_EDGES_SECTION) Node number out of range: %d",
                    j);
        if (i == j)
            Alg->eprintf("(FIXED_EDGES_SECTION) Illgal edge: %d to %d", i, j);
        Ni = &Alg->NodeSet[i];
        Nj = &Alg->NodeSet[Alg->ProblemType == LKH::ATSP ? j + Alg->Dimension / 2 : j];
        if (!Ni->FixedTo1)
            Ni->FixedTo1 = Nj;
        else if (!Ni->FixedTo2)
            Ni->FixedTo2 = Nj;
        else
            Alg->eprintf("(FIXED_EDGES_SECTION) Illegal fix: %d to %d", i, j);
        if (!Nj->FixedTo1)
            Nj->FixedTo1 = Ni;
        else if (!Nj->FixedTo2)
            Nj->FixedTo2 = Ni;
        else
            Alg->eprintf("(FIXED_EDGES_SECTION) Illegal fix: %d to %d", i, j);
        /* Cycle check */
        N = Ni;
        do {
            NNext = N->FixedTo1 != NPrev ? N->FixedTo1 : N->FixedTo2;
            NPrev = N;
            Count++;
        } while ((N = NNext) && N != Ni);
        if (N == Ni && Count != Alg->Dimension)
            Alg->eprintf("(FIXED_EDGES_SECTION) Illegal fix: %d to %d", i, j);
        if (!Alg->fscanint(Alg->ProblemFile, &i))
            i = -1;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
}

static void Read_NODE_COORD_SECTION(LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *N;
    int Id, i;

	CheckSpecificationPart(Alg);
    if (Alg->CoordType != LKH::TWOD_COORDS && Alg->CoordType != LKH::THREED_COORDS)
        Alg->eprintf("NODE_COORD_SECTION conflicts with NODE_COORD_TYPE: %s",
                Alg->NodeCoordType);
    if (!Alg->FirstNode)
		CreateNodes(Alg);
    N = Alg->FirstNode;
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    for (i = 1; i <= Alg->Dimension; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
            Alg->eprintf("Missing nodes in NODE_COORD_SECTION");
        if (Id <= 0 || Id > Alg->Dimension)
            Alg->eprintf("(NODE_COORD_SECTION) Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
            Alg->eprintf("(NODE_COORD_SECTION) Node number occours twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%lf", &N->X))
            Alg->eprintf("Missing X-coordinate in NODE_COORD_SECTION");
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Y))
           Alg-> eprintf("Missing Y-coordinate in NODE_COORD_SECTION");
        if (Alg->CoordType == LKH::THREED_COORDS
            && !fscanf(Alg->ProblemFile, "%lf", &N->Z))
            Alg->eprintf("Missing Z-coordinate in NODE_COORD_SECTION");
        if (Alg->Name && !strcmp(Alg->Name, "d657")) {
            N->X = (float) N->X;
            N->Y = (float) N->Y;
        }
    }
    N = Alg->FirstNode;
    do
        if (!N->V && N->Id <= Alg->Dimension)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
       Alg-> eprintf("(NODE_COORD_SECTION) No coordinates given for node %d",
                N->Id);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
}

static void Read_NODE_COORD_TYPE(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->NodeCoordType = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("NODE_COORD_TYPE: string expected");
    for (i = 0; i < strlen(Alg->NodeCoordType); i++)
        Alg->NodeCoordType[i] = (char) toupper(Alg->NodeCoordType[i]);
    if (!strcmp(Alg->NodeCoordType, "TWOD_COORDS"))
        Alg->CoordType = LKH::TWOD_COORDS;
    else if (!strcmp(Alg->NodeCoordType, "THREED_COORDS"))
        Alg->CoordType = LKH::THREED_COORDS;
    else if (!strcmp(Alg->NodeCoordType, "NO_COORDS"))
        Alg->CoordType = LKH::NO_COORDS;
    else
        Alg->eprintf("Unknown NODE_COORD_TYPE: %s", Alg->NodeCoordType);
}

static void Read_TOUR_SECTION(FILE ** File,LKH::LKHAlg *Alg)
{
    LKH::LKHAlg::Node *First = 0, *Last = 0, *N, *Na;
    int i, k;

    if (Alg->TraceLevel >= 1) {
        Alg->printff("Reading ");
        if (File == &Alg->InitialTourFile)
            Alg->printff("INITIAL_TOUR_FILE: \"%s\" ... ", Alg->InitialTourFileName);
        else if (File == &Alg->InputTourFile)
            Alg->printff("INPUT_TOUR_FILE: \"%s\" ... ", Alg->InputTourFileName);
        else if (File == &Alg->SubproblemTourFile)
            Alg->printff("SUBPROBLEM_TOUR_FILE: \"%s\" ... ",
                    Alg->SubproblemTourFileName);
        else
            for (i = 0; i < Alg->MergeTourFiles; i++)
                if (File == &Alg->MergeTourFile[i])
                    Alg->printff("MERGE_TOUR_FILE: \"%s\" ... ",
                            Alg->MergeTourFileName[i]);
    }
    if (!Alg->FirstNode)
		CreateNodes(Alg);
    N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (Alg->ProblemType == LKH::ATSP)
        Alg->Dimension /= 2;
    if (!Alg->fscanint(*File, &i))
        i = -1;
    for (k = 0; k <= Alg->Dimension && i != -1; k++) {
        if (i <= 0 || i > Alg->Dimension)
            Alg->eprintf("(TOUR_SECTION) Node number out of range: %d", i);
        N = &Alg->NodeSet[i];
        if (N->V == 1 && k != Alg->Dimension)
            Alg->eprintf("(TOUR_SECTION) Node number occours twice: %d", N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            if (Alg->ProblemType == LKH::ATSP) {
                Na = N + Alg->Dimension;
                Na->V = 1;
            } else
                Na = 0;
            if (File == &Alg->InitialTourFile) {
                if (!Na)
                    Last->InitialSuc = N;
                else {
                    Last->InitialSuc = Na;
                    Na->InitialSuc = N;
                }
            } else if (File == &Alg->InputTourFile) {
                if (!Na)
                    Last->InputSuc = N;
                else {
                    Last->InputSuc = Na;
                    Na->InputSuc = N;
                }
            } else if (File == &Alg->SubproblemTourFile) {
                if (!Na)
                    (Last->SubproblemSuc = N)->SubproblemPred = Last;
                else {
                    (Last->SubproblemSuc = Na)->SubproblemPred = Last;
                    (Na->SubproblemSuc = N)->SubproblemPred = Na;
                }
            } else {
                for (i = 0; i < Alg->MergeTourFiles; i++) {
                    if (File == &Alg->MergeTourFile[i]) {
                        if (!Na)
                            Last->MergeSuc[i] = N;
                        else {
                            Last->MergeSuc[i] = Na;
                            Na->MergeSuc[i] = N;
                        }
                    }
                }
            }
            Last = N;
        }
        if (k < Alg->Dimension)
            Alg->fscanint(*File, &i);
        if (k == Alg->Dimension - 1)
            i = First->Id;
    }
    N = Alg->FirstNode;
    do
        if (!N->V)
            Alg->eprintf("(TOUR_SECTION) Node is missing: %d", N->Id);
    while ((N = N->Suc) != Alg->FirstNode);
    if (File == &Alg->SubproblemTourFile) {
        do {
            if (N->FixedTo1 &&
                N->SubproblemPred != N->FixedTo1
                && N->SubproblemSuc != N->FixedTo1)
                Alg->eprintf("Fixed edge (%d, %d) "
                        "does not belong to subproblem tour", N->Id,
                        N->FixedTo1->Id);
            if (N->FixedTo2 && N->SubproblemPred != N->FixedTo2
                && N->SubproblemSuc != N->FixedTo2)
                Alg->eprintf("Fixed edge (%d, %d) "
                        "does not belong to subproblem tour", N->Id,
                        N->FixedTo2->Id);
        } while ((N = N->Suc) != Alg->FirstNode);
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
    if (Alg->ProblemType == LKH::ATSP)
        Alg->Dimension *= 2;
    if (Alg->TraceLevel >= 1)
        Alg->printff("done\n");
}

static void Read_TYPE(LKH::LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->Type = Copy(gStrtok_r(0, Delimiters,&Alg->savePtr))))
        Alg->eprintf("TYPE: string expected");
    for (i = 0; i < strlen(Alg->Type); i++)
        Alg->Type[i] = (char) toupper(Alg->Type[i]);
    if (!strcmp(Alg->Type, "TSP"))
        Alg->ProblemType = LKH::TSP;
    else if (!strcmp(Alg->Type, "ATSP"))
        Alg->ProblemType = LKH::ATSP;
    else if (!strcmp(Alg->Type, "SOP")) {
       Alg-> ProblemType = LKH::SOP;
        Alg->eprintf("TYPE: Type not implemented: %s", Alg->Type);
    } else if (!strcmp(Alg->Type, "HCP"))
        Alg->ProblemType = LKH::HCP;
    else if (!strcmp(Alg->Type, "CVRP")) {
        Alg->ProblemType = LKH::CVRP;
        Alg->eprintf("TYPE: Type not implemented: %s", Alg->Type);
    } else if (!strcmp(Alg->Type, "TOUR")) {
        Alg->ProblemType = LKH::TOUR;
        Alg->eprintf("TYPE: Type not implemented: %s", Alg->Type);
    } else if (!strcmp(Alg->Type, "HPP"))
        Alg->ProblemType = LKH::HPP;
    else
        Alg->eprintf("Unknown TYPE: %s", Alg->Type);
}

/*
   The ReadTour function reads a tour from a file.

   The format is as follows: 

   OPTIMUM = <real>
   Known optimal tour length. A run will be terminated as soon as a tour 
   length less than or equal to optimum is achieved.
   Default: MINUS_INFINITY.

   TOUR_SECTION :
   A tour is specified in this section. The tour is given by a list of integers
   giving the sequence in which the nodes are visited in the tour. The tour is
   terminated by a -1. 

   EOF
   Terminates the input data. The entry is optional.

   Other keywords in TSPLIB format may be included in the file, but they are 
   ignored.
*/

void LKH::LKHAlg::ReadTour(char *FileName, FILE ** File)
{
    char *Line, *Keyword, *Token;
    unsigned int i;
    int Done = 0;

    if (!(*File = fopen(FileName, "r")))
        eprintf("Cannot open tour file: \"%s\"", FileName);
    while ((Line = ReadLine(*File))) {
        if (!(Keyword = gStrtok_r(Line, Delimiters,&savePtr)))
            continue;
        for (i = 0; i < strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "COMMENT") ||
            !strcmp(Keyword, "DEMAND_SECTION") ||
            !strcmp(Keyword, "DEPOT_SECTION") ||
            !strcmp(Keyword, "DISPLAY_DATA_SECTION") ||
            !strcmp(Keyword, "DISPLAY_DATA_TYPE") ||
            !strcmp(Keyword, "EDGE_DATA_FORMAT") ||
            !strcmp(Keyword, "EDGE_DATA_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_FORMAT") ||
            !strcmp(Keyword, "EDGE_WEIGHT_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_TYPE") ||
            !strcmp(Keyword, "FIXED_EDGES_SECTION") ||
            !strcmp(Keyword, "NAME") ||
            !strcmp(Keyword, "NODE_COORD_SECTION") ||
            !strcmp(Keyword, "NODE_COORD_TYPE")
            || !strcmp(Keyword, "TYPE"));
        else if (strcmp(Keyword, "OPTIMUM") == 0) {
            if (!(Token = gStrtok_r(0, Delimiters,&savePtr)) ||
                !sscanf(Token, GainInputFormat, &Optimum))
                eprintf("[%s] (OPTIMUM): integer expected", FileName);
        } else if (strcmp(Keyword, "DIMENSION") == 0) {
            int Dim = 0;
            if (!(Token = gStrtok_r(0, Delimiters,&savePtr)) ||
                !sscanf(Token, "%d", &Dim))
                eprintf("[%s] (DIMENSION): integer expected", FileName);
            if (Dim != DimensionSaved)
                eprintf
                    ("[%s] (DIMENSION): does not match problem dimension",
                     FileName);
        } else if (!strcmp(Keyword, "TOUR_SECTION")) {
            Read_TOUR_SECTION(File,this);
            Done = 1;
        } else if (!strcmp(Keyword, "EOF"))
            break;
        else
            eprintf("[%s] Unknown Keyword: %s", FileName, Keyword);
    }
    if (!Done)
        eprintf("Missing TOUR_SECTION in tour file: \"%s\"", FileName);
    fclose(*File);
}
