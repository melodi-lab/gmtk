/*
 * GMTK_BoundaryTriangulate.h
 *   The GMTK Boundary Algorithm and Graph Triangulation Support Routines
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu> & Chris Bartels <bartels@ee.washington.edu>
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 * $Header$
 */

#ifndef GMTK_BOUNDARYTRIANGULATE_H
#define GMTK_BOUNDARYTRIANGULATE_H

#include <vector>
#include <string>
#include <list>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_RandomVariable.h"
#include "GMTK_FileParser.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_MaxClique.h"
#include "GMTK_Timer.h"

#include "debug.h"

// class mention for forward references.
class GraphicalModel;

class BoundaryTriangulate : public IM 
{
  friend class FileParser;
  friend class GraphicalModel;
  friend class GMTemplate;

private:

  // the file parser for this model.
  FileParser& fp;

  // number of chunks in which to find interface boundary
  const unsigned M;

  // chunk skip, Number of chunks that should exist between boundaries
  const unsigned S;

  // Structures over which the best boundary search is performed.
  struct BoundaryFindingStructures {
    // create sets P, C1, C2, C3, and E, from graph 
    // unrolled M+1 times .
    // prologue
    set<RandomVariable*> P;
    // 1st chunk, 1 chunk long
    set<RandomVariable*> C1;
    // 2nd chunk, M chunks long
    set<RandomVariable*> C2;
    // 3rd chunk, 1 chunk long
    set<RandomVariable*> C3;
    // epilogue
    set<RandomVariable*> E;
  };

  // Structures used to form partitions from a given
  // boundary. Once the partitions are formed, they 
  // can be triangulated.
  struct PartitionStructures {
    // create sets P, C1, C2, C3, and E, from graph 
    // unrolled M+S-1 times .
    // prologue
    set<RandomVariable*> P;
    // 1st chunk
    set<RandomVariable*> C1;
    // 2nd chunk
    set<RandomVariable*> C2;
    set<RandomVariable*> Cextra;
    // epilogue
    set<RandomVariable*> E;
  };

  ////////////////////////////////////////////////////////////
  // data types and such stuff
  ////////////////////////////////////////////////////////////

  enum BasicTriangulateHeuristic { /* S */ TH_MIN_SIZE = 1,        
				    /* T */ TH_MIN_TIMEFRAME = 2,
				    /* F */ TH_MIN_FILLIN = 3,
				    /* W */ TH_MIN_WEIGHT = 4,
				    // use average entropy in CPTs
				    /* E */ TH_MIN_ENTROPY = 5,
				    // use position in .str file
				    /* P */ TH_MIN_POSITION_IN_FILE = 6,
				    // use triangulation hint of variables
				    // in .str file given by user.
				    /* H */ TH_MIN_HINT = 7,
				    // use weight, but don't use any information
				    // about determinism of variables in clique.
				    /* N */ TH_MIN_WEIGHT_NO_D = 8,
				    // random triangulation
				    /* R */ TH_RANDOM = 9, 
				    // reverse time frame
				    /* X */ TH_MAX_TIMEFRAME = 10
  };
 
  typedef enum { 
    NO_EDGES,
    ALL_EDGES,
    RANDOM_EDGES,
    LOCALLY_OPTIMAL_EDGES
  } extraEdgeHeuristicType;

  enum TriangulateStyles { 
    TS_ANNEALING, 
    TS_EXHAUSTIVE,
    TS_MCS,
    TS_COMPLETED,
    TS_BASIC,
    TS_FRONTIER,
    TS_PRE_EDGE_ALL, 
    TS_PRE_EDGE_LO, 
    TS_PRE_EDGE_RANDOM, 
    TS_ELIMINATION_HEURISTICS, 
    TS_NON_ELIMINATION_HEURISTICS,
    TS_ALL_HEURISTICS 
  };

  struct TriangulateHeuristics {
    unsigned               numberTrials;
    TriangulateStyles      style;
    extraEdgeHeuristicType extraEdgeHeuristic;
    vector<BasicTriangulateHeuristic> heuristic_vector;
    unsigned numRandomTop;
    string basic_method_string;

    TriangulateHeuristics() {
      init();
    }
    void init() {
      // number of random re-tries
      numberTrials = 1;
      // default style
      style = TS_BASIC;
      // Add no extra edges 
      extraEdgeHeuristic = NO_EDGES;
      // default basic case:
      basic_method_string = "WFT";
      // first by weight
      heuristic_vector.push_back(TH_MIN_WEIGHT);
      // then by fill in if weight in tie
      heuristic_vector.push_back(TH_MIN_FILLIN);
      // and lastly by time frame (earliest first)
      heuristic_vector.push_back(TH_MIN_TIMEFRAME);
      // Default number of top nodes to randomly choose from
      // when eliminating. I.e., if we use heuristic 'W', we
      // take the top N of the nodes as ranked by 'W' and, rather
      // than just eliminating the top one, we randomly choose from
      // among the top N as the node to next eliminate.
      numRandomTop = 3;
    }
  };

  // TODO: change name to BoundaryHeuristic
  enum BoundaryHeuristic { /* S */ IH_MIN_SIZE = 1,        
			   /* F */ IH_MIN_FILLIN = 2,
			   /* W */ IH_MIN_WEIGHT = 3,
			   /* N */ IH_MIN_WEIGHT_NO_D = 4,
			   /* E */ IH_MIN_ENTROPY = 5,
			   /* M */ IH_MIN_MAX_CLIQUE = 6,
			   /* C */ IH_MIN_MAX_C_CLIQUE = 7,
			   /* A */ IH_MIN_STATE_SPACE = 8,
			   /* Q */ IH_MIN_C_STATE_SPACE = 9,
			   /* p */ IH_MIN_MIN_POSITION_IN_FILE = 10,
			   /* P */ IH_MIN_MAX_POSITION_IN_FILE = 11,
			   /* t */ IH_MIN_MIN_TIMEFRAME = 12,
			   /* T */ IH_MIN_MAX_TIMEFRAME = 13,
			   /* h */ IH_MIN_MIN_HINT = 14,
			   /* H */ IH_MIN_MAX_HINT = 15
  };

  // While much of the code for finding the best right and best left
  // interface boundaries is symmetric, there are a few places where
  // the code needs to do something different. Rather than pass in an
  // argument deep inside a set of routine calls, we keep this member
  // variable, set in the outer loops, to say if we are currently
  // finding the left interface boundary (set to true) or finding the right
  // interface boundary (set to false).
  bool findingLeftInterface;

  // Support for unTriangulating partitions. This stuff
  // could go into the Partition class, but we don't want
  // to include all the STL code for triangulation/untriangulting
  // since Partitiona and GMTemplate will most often be used
  // for inference, and where a triangulation will simply
  // come from a set of maxcliques.
  typedef pair<RandomVariable*, set<RandomVariable*> > nghbrPairType; 

  // original neighbors of nodes so that deTriangulate() will work.
  vector<nghbrPairType>  orgnl_P_nghbrs;
  vector<nghbrPairType>  orgnl_C_nghbrs;
  vector<nghbrPairType>  orgnl_E_nghbrs;

  // The timer for the anytime algorithm. If this variable
  // is non-NULL, some of the routines will check 'timer' and if
  // it has expired, the routines will return. If this variable is NULL,
  // there will be no effect.
  TimerClass* timer;

  // variable to keep track of recursion depth of boundary algorithm,
  // for debugging and informational purposes.
  unsigned boundaryRecursionDepth;


  /////////////////////////////////////////////////////
  // Private support routines
  /////////////////////////////////////////////////////

  // delete a set of variables
  void deleteNodes(const set<RandomVariable*>& nodes);

  // given a string, create a vector of triangulation heuristics
  void parseTriHeuristicString(const string& tri_heur_str,TriangulateHeuristics& th);
  void createVectorBoundaryHeuristic(const string& th,
				     vector<BoundaryHeuristic>& th_v);

  // computes the fill in of a set of variables.
  int computeFillIn(const set<RandomVariable*>& nodes);

  // compute the weight of a vector of cliques
  double graphWeight(vector<MaxClique>& cliques);
  double graphWeight(vector<MaxClique>& cliques, 
		     // true if we should use JT rather than normal weight.
		     const bool useJTWeight,
		     // if useJTWeight is true, this gives nodes that
		     // root must cover.
		     const set<RandomVariable*>& interfaceNodes);
  
  void fillAccordingToCliques(
    const vector<MaxClique>& cliques 
  );


  ////////////////////////////////////////////////////////////
  // options
  ////////////////////////////////////////////////////////////
  // don't memoize the boundary to save memory.
  bool noBoundaryMemoize;
  // fraction of boundary to traverse (randomly), produces sub-optimal boundary but
  // runs much faster.
  double boundaryTraverseFraction;

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // The Main Routines To Form optimal P,C,E partitions
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  void findInterfacePartitions(
   // input params
   const set<RandomVariable*>& P_u1,
   const set<RandomVariable*>& C1_u1,
   const set<RandomVariable*>& Cextra_u1, // non-empty only when S > M
   const set<RandomVariable*>& C2_u1,
   const set<RandomVariable*>& E_u1,
   map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
   map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
   const set<RandomVariable*>& left_C_l_u2C2,
   const set<RandomVariable*>& C_l_u2C2,
   // output params
   GMTemplate& gm_template);



  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // The Main Triangulation Routines
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  // when triangulating a partition, jt-weight likes
  // to know the nodes in the left and/or right partition (if any
  // are availalbe). These variables, if not set to NULL, contain
  // pointers to these. These are set in the main public triangulation routine.
  // Left-partition nodes (NULL if no left partition relative to current partition)
  set<RandomVariable*> *lp_nodes;
  // right partition nodes.
  set<RandomVariable*> *rp_nodes;  

  // Calls method which triangulates once, support routine for triangulate 
  void triangulateOnce(// input: nodes to be triangulated
  		       const set<RandomVariable*>& nodes,
                       // use JT weight rather than sum of weight
		       const bool jtWeight,
		       // nodes that a JT root must contain (ok to be empty).
		       const set<RandomVariable*>& nodesRootMustContain,
		       // triangulation heuristic method
		       const TriangulateHeuristics& tri_heur,
		       // original neighbor structures
		       vector<nghbrPairType>& orgnl_nghbrs,
		       // output: resulting max cliques
		       vector<MaxClique>& best_cliques,
		       // output: string giving resulting method used
		       string& meth_str );
  
  // High-level generic graph triangulation using optionally all methods below.
  void triangulatePartition(// input: nodes to be triangulated
			    const set<RandomVariable*>& nodes,
			    // use JT weight rather than sum of weight
			    const bool jtWeight,
			    // nodes that a JT root must contain (ok to be empty).
			    const set<RandomVariable*>& nodesRootMustContain,
			    // triangulation heuristic method
			    const TriangulateHeuristics& tri_heur,
			    // original neighbor structures
			    vector<nghbrPairType>& orgnl_nghbrs,
			    // output: resulting max cliques
			    vector<MaxClique>& best_cliques,
			    // output: string giving resulting method used
			    string& best_meth_str,
			    // weight to best
			    double& best_weight,
			    // Prefix for new best_meth_str 
			    string  best_method_prefix = "");


  // version of triangulatePartition() above that takes triangulation
  // heuristic strings.
  void triangulatePartition(// input: nodes to be triangulated
			    const set<RandomVariable*>& nodes,
			    // use JT weight rather than sum of weight
			    const bool jtWeight,
			    // nodes that a JT root must contain (ok to be empty).
			    const set<RandomVariable*>& nodesRootMustContain,
			    // triangulation heuristic method
			    const string& tri_heur_str,
			    // original neighbor structures
			    vector<nghbrPairType>& orgnl_nghbrs,
			    // output: resulting max cliques
			    vector<MaxClique>& best_cliques,
			    // output: string giving resulting method used
			    string& best_meth_str,
			    // weight to best
			    double& best_weight, 
			    // Prefix for new best_meth_str 
			    string  best_method_prefix = "" ) 
  {
    TriangulateHeuristics tri_heur;
    parseTriHeuristicString(tri_heur_str,tri_heur);
    triangulatePartition(nodes,jtWeight,nodesRootMustContain,
			 tri_heur,orgnl_nghbrs,
			 best_cliques,best_meth_str,best_weight,best_method_prefix);
  }


  // include version that always uses the standard weight measure
  // (i.e., does not use jt-weight).
  void triangulatePartitionWeight(// input: nodes to be triangulated
				  const set<RandomVariable*>& nodes,
				  // triangulation heuristic method
				  const TriangulateHeuristics& tri_heur,
				  // original neighbor structures
				  vector<nghbrPairType>& orgnl_nghbrs,
				  // output: resulting max cliques
				  vector<MaxClique>& best_cliques,
				  // output: string giving resulting method used
				  string& best_meth_str,
				  // weight to best
				  double& best_weight)
  {
    const set <RandomVariable*> emptySet;
    triangulatePartition(nodes, false, emptySet, tri_heur, orgnl_nghbrs,
			 best_cliques, best_meth_str, best_weight);
  }


  // The basic triangulation heuristic routine. Given a set of nodes
  // (that have valid 'neighbors' members, triangulate it using the
  // heuristic(s) given and return the resuting set of cliques in no
  // particular order (i.e., not in RIP order).  Note that more than
  // one heuristic can be used a a time in order to break ties. As a
  // last resort (if all the heuristics agree on the same next set of
  // nodes to eliminate) then a uniformly at random choice is
  // made. This routine is guaranteed to run fast an so can be used in
  // online mode. Also, this routine can be called multiple times
  // where the caller chooses the results with the best cliques --
  // this is because each call might produce a different clique set
  // via the internal randomness that can occur if a tie occurs.
  void basicTriangulate(const set<RandomVariable*>& nodes,
			const vector<BasicTriangulateHeuristic>& th_v,
			const unsigned numRandomTop,
			vector<RandomVariable*>& orderedNodes,
			vector<MaxClique>& cliques,
			const bool findCliques = true);

  // triangulate by simulated annealing
  void triangulateSimulatedAnnealing(
    const set<RandomVariable*>& nodes,
    const bool                  jtWeight,
    const set<RandomVariable*>& nodesRootMustContain,
    vector<MaxClique>&          best_cliques,
    vector<RandomVariable*>&    best_order,
    string&                     comment 
    );

  // Triangulate by maximum cardinality search
  void triangulateMaximumCardinalitySearch( 
    const set<RandomVariable*>& nodes,
    vector<MaxClique>&          cliques,
    vector<RandomVariable*>&    order
    );

  // This procedure is like triangulateMaximumCardinalitySearch except
  // that it tests if the original graph was triangulated
  bool triangulateMCSIfNotTriangulated( 
    const set<RandomVariable*>& nodes,
    vector<MaxClique>&          cliques
  );

  // Check if graph is chordal 
  bool chordalityTest( 
    const set<RandomVariable*>& nodes
  );

  // Check chordality and get cliques in RIP order
  bool getCliques( 
    const set<RandomVariable*>& nodes,
    vector<MaxClique>&          cliques
  );

  // triangulation by simple completion
  void triangulateCompletePartition(const set<RandomVariable*>& nodes,
				    vector<MaxClique>&          cliques
				    );

  // triangulation by frontier algorithm
  void triangulateFrontier(const set<RandomVariable*>& nodes,
			   vector<MaxClique>&          cliques
			   );


  // triangulate by exhaustive search, takes a *LONG* time.
  void triangulateExhaustiveSearch(const set<RandomVariable*>&  nodes,
				   const bool jtWeight,
				   const set<RandomVariable*>& nodesRootMustContain,
				   const vector<nghbrPairType>& orgnl_nghbrs,
				   vector<MaxClique>&           cliques
				   );

  // Triangulate by pre-specified elimination order
  void triangulateElimination(// input: nodes to be triangulated
			      const set<RandomVariable*> nodes,
			      // elimination order 
			      vector<RandomVariable*> orderedNodes,  
			      // output: resulting max cliques
			      vector<MaxClique>& cliques);


  ////////////////////////////////////////////////////////////
  // triangulate using elimination with a number of basic heuristics, 
  // returning the best.
  double tryEliminationHeuristics(
    const set<RandomVariable*>& nodes,
    const bool                  jtWeight,
    const set<RandomVariable*>& nodesRootMustContain,
    vector<nghbrPairType>&      orgnl_nghbrs,
    vector<MaxClique>&          best_cliques,
    string&                     best_method,
    double&                     best_weight,
    string                      best_method_prefix = ""
    );

  ////////////////////////////////////////////////////////////
  // triangulate using using a number of non-elimination based heuristics, 
  // returning the best.
  double tryNonEliminationHeuristics(
    const set<RandomVariable*>& nodes,
    const bool                  jtWeight,
    const set<RandomVariable*>& nrmc,         // nrmc = nodes root must contain
    vector<nghbrPairType>&      orgnl_nghbrs,
    vector<MaxClique>&          best_cliques,
    string&                     best_method,
    double&                     best_weight
    );

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // Boundary Algorithm Routines
  ////////////////////////////////////////////////////////////

  // general routine to compute the score of a candidate interface.
  void interfaceScore(
	 const vector<BoundaryHeuristic>& bnd_heur_v,
	 const set<RandomVariable*>& C_l,
	 // more input variables for use when MIN_CLIQUE is active
	 const set<RandomVariable*>& left_C_l,
	 const TriangulateHeuristics& tri_heur,
	 const set<RandomVariable*>& P_u1,
	 const set<RandomVariable*>& C1_u1,
	 const set<RandomVariable*>& Cextra_u1,
	 const set<RandomVariable*>& C2_u1,
	 const set<RandomVariable*>& E_u1,
	 // these next 2 should be const, but there is no "op[] const"
	 map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
	 map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1,
	 // output score
	 vector<float>& score);

  // exponential time routines to find best left/right interfaces
  // and corresponding partitions. Since the operations
  // for finding the best left and right interfaces are symmetric,
  // the determiniation of if we are searching for the best
  // left interface or right interface is determined entirely
  // based on the arguments that are passed into these routines. 
  // Note that for simplicity, the names have been defined
  // in terms of the left interface, but that is only for ease of
  // understanding.
  void findBestInterface(
	     const set<RandomVariable*> &C1,
	     const set<RandomVariable*> &C2,
	     const set<RandomVariable*> &C2_1,
	     const set<RandomVariable*> &C3,
	     set<RandomVariable*> &left_C_l,
	     set<RandomVariable*> &C_l,
	     vector<float>& best_score,
	     const vector<BoundaryHeuristic>& bnd_heur_v,
	     const bool recurse,
	     const TriangulateHeuristics& tri_heur,
	     const set<RandomVariable*>& P_u1,
	     const set<RandomVariable*>& C1_u1,
	     const set<RandomVariable*>& Cextra_u1,
	     const set<RandomVariable*>& C2_u1,
	     const set<RandomVariable*>& E_u1,
	     // these next 2 should be const, but there is no "op[] const"
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1);

  void findBestInterface(
             const set<RandomVariable*> &left_C_l,
	     const set<RandomVariable*> &C_l,
	     const set<RandomVariable*> &C2,
	     const set<RandomVariable*> &C2_1,
	     const set<RandomVariable*> &C3,
	     set< set<RandomVariable*> >& setset,
	     set<RandomVariable*> &best_left_C_l,
	     set<RandomVariable*> &best_C_l,
	     vector<float>& best_score,
	     const vector<BoundaryHeuristic>& bnd_heur_v,
	     const TriangulateHeuristics& tri_heur,
	     const set<RandomVariable*>& P_u1,
	     const set<RandomVariable*>& C1_u1,
	     const set<RandomVariable*>& Cextra_u1,
	     const set<RandomVariable*>& C2_u1,
	     const set<RandomVariable*>& E_u1,
	     // these next 2 should be const, but there is no "op[] const"
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C1_u1,
	     map < RandomVariable*, RandomVariable* >& C2_u2_to_C2_u1);


  //////////////////////////////////////////////////////////////////////////// 
  // Custom data structures for fast Maximum Cardinality Search, fill-in 
  // computation, and chordality test. 
  //////////////////////////////////////////////////////////////////////////// 
  class triangulateNode;
  class triangulateNodeList;

  typedef vector<triangulateNode*> triangulateNeighborType;

  class triangulateNode {
 
    friend class triangulateNodeList;

    public:
      triangulateNode();
      triangulateNode(RandomVariable* random_variable);

      RandomVariable*          randomVariable;
      triangulateNeighborType  neighbors;
      triangulateNodeList*     nodeList; 
      unsigned                 cardinality;
      unsigned                 position;
      bool                     eliminated;
      bool                     marked;
      triangulateNeighborType  parents;
      triangulateNeighborType  nonChildren;
    
    private:
      triangulateNode*         previousNode;
      triangulateNode*         nextNode;
  };

  class triangulateNodeList {

    public:

      triangulateNodeList(); 
      void             push_back(triangulateNode* node); 
      triangulateNode* pop_back(); 
      void             erase(triangulateNode* node);
      unsigned size() { return(list_length); }
      triangulateNode* back() { return(last); }

      triangulateNode* operator[] (unsigned i);
  
    private:

      triangulateNode* last; 
      unsigned         list_length;
  }; 

  void fillTriangulateNodeStructures( 
    const set<RandomVariable*>& orgnl_nodes,
    vector<triangulateNode>&    new_nodes 
  );

  //////////////////////////////////////////////////////////////////////////// 
  // Edge class
  //////////////////////////////////////////////////////////////////////////// 
  class edge {
    private:
      triangulateNode* first_node;
      triangulateNode* second_node;

      void assign( triangulateNode* a, triangulateNode* b) {
        if (a < b) {
          first_node  = a;
          second_node = b;
        } 
        else {
          first_node  = b;
          second_node = a;
        } 
      }

    public:
  
      triangulateNode* first() const { return(first_node); };
      triangulateNode* second() const { return(second_node); };

      edge() { first_node = NULL; second_node = NULL; }

      edge operator= (const edge& e) {
        assign( e.first(), e.second() ); 
        return(*this);
      }

      edge( triangulateNode* a, triangulateNode* b) {
        assign(a,b); 
      }

      bool operator< (const edge& e) const {
        if (first_node != e.first()) {
          return(first_node<e.first());
        }
        else { 
          return(second_node<e.second());
        }
      }
  };
 
  //////////////////////////////////////////////////////////////////////////// 
  // O(n+e) routines for Maximum Cardinality Search, fill-in 
  // computation, and chordality test. 
  //////////////////////////////////////////////////////////////////////////// 
  void maximumCardinalitySearch(
    vector<triangulateNode>&         nodes,
    list<vector<triangulateNode*> >& cliques,
    vector<triangulateNode*>&        order,
    bool                             randomize_order 
    );

  void fillInComputation(
    vector<triangulateNode*>& ordered_nodes 
  );
 
  bool testZeroFillIn( 
    vector<triangulateNode*>& ordered_nodes 
  );

  void listVectorCliquetoVectorSetClique(
    const list<vector<triangulateNode*> >& lv_cliques,
    vector<MaxClique>&                     vs_cliques
  );

  //////////////////////////////////////////////////////////////////////////// 
  // Overloaded functions for saving and restoring graph structure when using 
  // triangulateNodes
  //////////////////////////////////////////////////////////////////////////// 

  void saveCurrentNeighbors(
    const set<RandomVariable*> nodes,
    vector<nghbrPairType>& orgnl_nghbrs
  );

  void saveCurrentNeighbors(
    Partition &prt,
    vector<nghbrPairType>& orgnl_nghbrs) 
  {
    saveCurrentNeighbors(prt.nodes,orgnl_nghbrs);
  }

  void restoreNeighbors(vector<nghbrPairType>& orgnl_nghbrs);

  typedef pair<triangulateNode*, vector<triangulateNode*> > 
    triangulateNghbrPairType; 

  void saveCurrentNeighbors(
    vector<triangulateNode>&          nodes,
    vector<triangulateNghbrPairType>& orgnl_nghbrs
    );

  void restoreNeighbors(
    vector<triangulateNghbrPairType>& orgnl_nghbrs
  );

  //////////////////////////////////////////////////////////////////////////// 
  // Extra edges 
  //////////////////////////////////////////////////////////////////////////// 
 
  void fillParentChildLists(
    vector<triangulateNode>& nodes
  );

  void addExtraEdgesToGraph(
    const set<RandomVariable*>&  nodes,
    const extraEdgeHeuristicType edge_heuristic
  );

  void addExtraEdgesToGraph(
    vector<triangulateNode>&     nodes,
    const extraEdgeHeuristicType edge_heuristic
  );

  void addEdgesToNode(
    triangulateNeighborType& parent_set,  
    triangulateNode* const       child, 
    triangulateNode* const       grandchild, 
    const extraEdgeHeuristicType edge_heuristic,
    vector<edge>&                extra_edges 
  );

  void addEdges(
    const vector<edge>& extra_edges
  );

  //////////////////////////////////////////////////////////////////////////// 
  // Support routine for simulated annealing
  //////////////////////////////////////////////////////////////////////////// 
  unsigned annealChain(
    vector<triangulateNode>&  nodes,
    const bool jtWeight,
    const set<RandomVariable*>& nodesRootMustContain,
    vector<triangulateNode*>& crrnt_order,
    vector<triangulateNode*>& triangulate_best_order,
    double&                   best_graph_weight,
    double&                   best_this_weight,
    double                    temperature,
    unsigned                  iterations,
    double&                   weight_sum,         
    double&                   weight_sqr_sum,         
    vector<triangulateNghbrPairType>&    orgnl_nghbrs
    );

public:

  // Public Interface

  ////////////////////////////////////////////////////////////
  // constructors/destructors and other misc.
  ////////////////////////////////////////////////////////////
  BoundaryTriangulate(FileParser& arg_fp, 
		      const unsigned arg_M,
		      const unsigned arg_S,
		      double arg_boundaryTraverseFraction = 1.0) 
    : fp(arg_fp),M(arg_M),S(arg_S),
      noBoundaryMemoize(false),
      boundaryTraverseFraction(arg_boundaryTraverseFraction)
  {   
    assert ( M >= 1 ); assert( S >= 1); 
    timer = NULL;  // disable anytime by setting timer to NULL.
  }
  ~BoundaryTriangulate() {}

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  // Main interface Support for graph Triangulation.
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  // Main external interface to graph triangulation routine using a
  // simple textual string to determine the heuristics. This routine
  // will triangulate an entire template at once, putting cliques in
  // the template, but it includes options to only triangulate one
  // partition at a time if so desired.
  void triangulate(const string& tri_heur_str,
		   const bool jtWeight,
		   GMTemplate& gm_template,
		   bool doP = true,  // triangulate P
		   bool doC = true,  // triangulate C
		   bool doE = true   // triangulate E
		   );

  // A simple one-stop shop for good anytime algorithm triangulation,
  // runs only for a given amount of time.
  void anyTimeTriangulate(GMTemplate& gm_template,
			  const bool jtWeight,
			  bool doP = true, bool doC = true, bool doE = true);

  // Given the template, just unroll it flat-out a given number of
  // times and triangulate the result (possibly unconstrained
  // but depending on the heuristics given in 'th').
  void unrollAndTriangulate(const string& tri_heur_str,
			    const unsigned numTimes);


  // Given the template, compute the best partitions
  // using the heuristics that are provided. If 'findBestBoundary'
  // is true, this will run the exponential algorithm (which
  // isn't necessarily bad since for many structures
  // it is still efficient, but for some it might blow
  // up in amount of time to complete, so caller should
  // be forwarned)
  void findPartitions(const string& boundaryHeuristic,
		      const string& forceLeftRight,
		      const string& triHeuristic,
		      const bool findBestBoundary,
		      GMTemplate& gm_template);


  // use the timer given by arg, returning the old timer.
  // If no argument given, then don't use any timer.
  TimerClass* useTimer(TimerClass* arg = NULL) { 
    TimerClass* tmp = timer;
    timer = arg;
    return tmp;
  }

  // Ensure that the partitions in the given template are chordal, and
  // die with an error if not.
  void ensurePartitionsAreChordal(GMTemplate& gm_template);

  // interface to private boolean.
  void dontMemoizeBoundary() {
    noBoundaryMemoize = true;
  }

};


#endif

