/*-
 * GMTK_LatticeADT.h
 *      .h file for GMTK_LatticeADT.cc, HTK lattice support
 *      distributions.
 *
 *  Written by Gang Ji <gang@ee.washington.edu> and Jeff Bilmes <bilmes@uw.edu>
 * 
 *  $Header$
 * 
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */


#ifndef GMTK_LATTICE_ADT_H
#define GMTK_LATTICE_ADT_H


#include "GMTK_NamedObject.h"
#include "GMTK_Vocab.h"
#include "fileParser.h"
#include "shash_map_iter.h"
#include "logp.h"


/**
 * HTK lattice support
 */
class LatticeADT : public NamedObject {
 public:

  // set to true if the lattice nodes cpt for (node,node) use a score
  // which is the max of the score over all edges for a given
  // (node,node), or if it should return one. In the case of max, then
  // the lattice edges must divide out this max, so this variable
  // is used by both the lattice edge and lattice node cpts.
  static bool _latticeNodeUseMaxScore;

  /** frame rate **/
  static double _defaultFrameRate;

  LatticeADT();
  ~LatticeADT();

  // read from HTK lattice format file
  void readFromHTKLattice(iDataStreamFile &ifs, const Vocab &vocab);

  // read from GMTK master file
  void read(iDataStreamFile &is);
  
  // read the score options and return them encoded in an int.
  unsigned readScoreOptions(iDataStreamFile &is);
  void printScoreOptions(FILE* f);

  // this lattice is iterable
  inline bool iterable() const { return _latticeFile != NULL; }

  void seek(unsigned nmbr);
  void initializeIterableLattice(const string &fileName);
  void beginIterableLattice();
  void nextIterableLattice();

  // reset the frame indices
  void resetFrameIndices(unsigned numFrames);
  void setGMTKScores();
  void printLatticeInfo(FILE*f);

  inline bool useTimeParent() const {
    return _timeCardinality != 0;
  }

  friend class LatticeNodeCPT;
  friend class LatticeEdgeCPT;

 protected:
  /**
   * lattice edge
   */
  struct LatticeEdge {
    /** emission id from current state to next state */
    unsigned emissionId;
    /** acoustic score */
    logpr ac_score;
    /** acoustic score */
    logpr lm_score;
    /** posterior */
    logpr posterior;
    /** score used in GMTK */
    logpr gmtk_score;

    LatticeEdge() : emissionId(0), gmtk_score(1.0) {}
  };

  /*
   * For every pair of connected, nodes there can be one or more
   * edges.  This list contains those edges associated with a
   * connected pair of nodes. This also contains the precomputed score
   * to be used for the node->node transition for this set of edges
   * (which can either be unity, or can be the max score of each of
   * the edges) -- depending on this score, the lattice node will need
   * to be scored differently.
   */
  struct LatticeEdgeList {
    unsigned num_edges;
    // use an sArray without a destructor, so that when hash table
    // resizes, it doesn't delete the memory used.
    sArray_nd < LatticeADT::LatticeEdge > edge_array;

    /** max score used in GMTK for all the edges */
    logpr max_gmtk_score;

    LatticeEdgeList() : num_edges(0) {} 
  };

  /**
   * lattice node information
   */
  struct LatticeNode {
    /** absolute time for this node */
    float time;
    /** starting frame number */
    unsigned startFrame;
    /** ending frame number */
    unsigned endFrame;
    /** possible out-going edges */
    shash_map_iter<unsigned, LatticeEdgeList > edges;

    LatticeNode() : startFrame(0), endFrame(0), edges(shash_map_iter<unsigned, LatticeEdgeList>(1)) {}
    ~LatticeNode(); 
  };


  /**
   * renormalize posterior
   */
  void normalizePosterior();

  /** lattice nodes */
  LatticeNode *_latticeNodes;
  /** number of nodes in lattice this should be smaller than node cardinality */
  unsigned _numberOfNodes;
  /** number of links in lattice */
  unsigned _numberOfLinks;
  /** start node id */
  unsigned _start;
  /** end node id */
  unsigned _end;
  /** language model scale */
  double _lmscale;
  /** word penalty */
  double _wdpenalty;
  /** acoustic model scale */
  double _acscale;
  /** log base **/
  double _base;


  // frameRate is the frame rate (e.g., in units of frames per second,
  // if the time marks in the lattice are in seconds).
  double _frameRate;
  /*
   * How many frames of "relaxation" is allowed in CPT, both to the
   * left and to the right of the lattice node transition.
   */
  unsigned _frameRelax;

  /** if this is an iterable cpt */
  iDataStreamFile* _latticeFile; // the file pointer
  string _latticeFileName; // the file name
  unsigned _numLattices; // number of lattices in this file
  int _curNum; // the current lattice number
  string _curName; // the current lattice cpt name

  // the following is GM paramters
  unsigned _nodeCardinality;
  unsigned _wordCardinality;
	
  // lattice can optionally have time as parent.  In this
  // case, we need time cardinality. 0 mean no time parent.
  unsigned _timeCardinality;

  unsigned score_options;
};


#endif

