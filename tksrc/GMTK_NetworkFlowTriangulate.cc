/*
 * GMTK_NetworkFlowTriangulate.cc
 *   The GMTK Network flow triangulation support routines.
 *
 * Written by Mukund Narasimhan <mukundn@ee.washington.edu> 
 *
 * Copyright (c) 2005, < fill in later >
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

#include "debug.h"
#include "GMTK_NetworkFlowTriangulate.h"


void GMTK2Network::init(float (*computeWeight)(const set<RV*>&), 
			const std::set<RV*>& cLeft, 
			const std::set<RV*>& cMiddle, 
			const std::set<RV*>& cRight) {

  // 2 nodes for each RV in cMiddle + source node + terminal node
  int numNodes = cMiddle.size() * 2 + 2;      
  Network::init(numNodes);
  setSource(0); setTerminal(1);
  

  // Nodes 0 is the source, Node 1 is the terminal. 
  // Odd numbered vertices represent the "in" vertex and
  // Even numbered vertices represent the "out" vertex.
  // For each random variable, there is a single edge going 
  // from the "in" vertex to the "out" vertex and the capacity
  // of this edge is the weight of the random variable.
  int p=2; std::set<RV*>::iterator i,begin,end;


  begin = cMiddle.begin(); end = cMiddle.end();  
  for(i=begin;i!=end;++i) {
    _translation[*i] = p;
    _reverseTranslation[p] = *i;
    _reverseTranslation[p+1] = *i;
    p += 2;
  }
  
  // Set up edges from source to rv's in left interface.
  begin=cLeft.begin(); end=cLeft.end();
  for(i=begin;i!=end;++i) {
    // sanity check.
    std::set<RV*>::iterator k = cMiddle.find(*i);
    if (k!=end) {
      addDirectedEdge(0, inVertex(*i), -1);
    }
  }
  
  // Set up edges from rv's in right interface to terminal.
  begin=cRight.begin();end=cRight.end();
  for(i=begin;i!=end;++i) {
    // sanity check.
    std::set<RV*>::iterator k = cMiddle.find(*i);
    if (k!=end) {
      addDirectedEdge(outVertex(*i), 1, -1);
    }
  }
  
  // Set up remaining edges.
  begin = cMiddle.begin(); end = cMiddle.end();  
  for(i=begin;i!=end;++i) {
    RV* currRV = *i;
    std::set<RV*> rv; rv.insert(*i);
    float weight = computeWeight(rv);
    addDirectedEdge(inVertex(*i), outVertex(*i), weight);

    std::set<RV*>::iterator j,neighboursBegin,neighboursEnd;
    neighboursBegin = currRV->neighbors.begin();
    neighboursEnd = currRV->neighbors.end();
    for(j=neighboursBegin;j!=neighboursEnd;++j) {
      // Add an edge to only those neighbours which are present in cMiddle
      // if (i,j) is an edge in the DBN, then there will be an edge from
      // out(i) to in(j) in this network.
      std::set<RV*>::iterator k = cMiddle.find(*j);
      if (k!=end) {
	const int f = outVertex(*i); const int t = inVertex(*j);
	addDirectedEdge(f, t, -1);
      }
    }
  }

}




void GMTK2Network::findBoundary(std::set<RV*>& boundary, std::set<RV*>& preBoundary, bool forward) {
  EdgeSet::iterator ei,cutEnd;
  infoMsg(IM::Info, "Computing the boundary");
  ei=_minCut.begin(); cutEnd=_minCut.end();
  for(;ei!=cutEnd;++ei) {
    int node= (*ei)->sourceNode();
    if (node != 0) {
      boundary.insert(_reverseTranslation[(*ei)->sourceNode()]);
    }
  }

  infoMsg(IM::Info, "Computing the frontier.");
  NodeSet::iterator ni,reachableEnd;
  reachableEnd=_reachableSet.end();
  std::set<RV*>::iterator boundaryEnd = boundary.end();
  for(ni=_reachableSet.begin();ni!=reachableEnd;++ni) {
    int currNode = *ni;
    if (currNode != 0) {
      RV* rv = _reverseTranslation[currNode];
      if (boundary.find(rv) == boundaryEnd)
	preBoundary.insert(rv);
    }
  }
}


