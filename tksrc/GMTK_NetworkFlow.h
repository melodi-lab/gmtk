/*
 * GMTK_NetworkFlow.h
 *   The GMTK Network flow support routines.
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

#ifndef GMTK_NETWORKFLOW_H
#define GMTK_NETWORKFLOW_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <numeric>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "debug.h"

namespace networkFlow {
  extern const bool DebugMode;
  ////////////////////////////////////////////////////////////////////
  // The NetworkFlowException class is used for exception handling. //
  // All methods in this file only throw exceptions of type         //
  // NetworkFlowException                                           //
  ////////////////////////////////////////////////////////////////////
  class NetworkFlowException {    
  public:
    std::string  _errorString;
    int          _errNum;


    // Standard Error Types
    enum ERRORS { INTERNAL, 
		  FLOW_EXCEEDS_CAPACITY, 
		  NEGATIVE_FLOW, 
		  UNBOUNDED_FLOW,
		  INVALID_PATH,
		  EDGE_ALREADY_EXISTS,
		  BIDIRECTIONAL_FLOW,
		  EMPTY_PATH,
                  UNSUPPORTED_FORMAT,
                  FILE_NOT_FOUND
                  };

    

    NetworkFlowException(std::string s, int errNum) { 
      _errorString = s;  _errNum = errNum;
    }


    NetworkFlowException(bool internal, int errNum) {
      _errorString = errorMessage(internal, errNum); _errNum = errNum;
      std::cout << "Throwing exception : " << _errorString << std::endl;
    }

    
    std::string errorMessage(bool internal, int errNum);
  };








  ////////////////////////////////////////////////////////////////////
  // The NetworkEdge<NodeType> class encapsulates flow related      //
  // information for directed edges including capacities, residual  //
  // capacities, flows and costs.                                   //
  // It is expected that NodeType will be a lightweight object, like//
  // and int or a pointer.                                          //
  // The edge is inherently directed. Undirected edges are to be    //
  // simulated using two directional edges.                         // 
  ////////////////////////////////////////////////////////////////////
  template<class NodeType>
    class NetworkEdge {
    public:
    friend class NetworkPath;
    // Argument Summary.
    // ns       : source node
    // nt       : target node
    // cap      : capacity of the edge. A negative capacity indicates 
    //            infinite or unbounded capacity from ns to nt.
    // cst      : cost of the edge (currently unused).
    // crntFlow : current (or pre-existing) flow.
    NetworkEdge(NodeType ns, NodeType nt,
		double cap = -1, // std::numeric_limits<double>::max(), 
		double cst = 0, double crntFlow = 0) :
      _nodeS(ns), _nodeT(nt) {
      _capacity = cap;
      _cost = cst;
      _flow = crntFlow;
    }



    double capacity(void) const { return _capacity; }

    double cost(void) const { return _cost; }

    double flow(void) const { return _flow; }

    bool hasUnlimitedCapacity(void) const { return (_capacity < 0); }
    
    NodeType sourceNode(void) const { return _nodeS; }

    NodeType targetNode(void) const { return _nodeT; }

    double residualCapacity(void) const;

    bool capacityExists(void) const;
    
    void addFlow(double incFlow) ;

    void clearFlow(void);

    protected:
    NodeType               _nodeS;         // Source Node
    NodeType               _nodeT;         // Terminal Node 
    double                 _capacity;      // Edge Capacity. Negative value 
                                           // indicates infinite capacity.
    double                 _flow;          // Current flow on the edge.
    double                 _cost;          // Cost of the edge (currently unused).



    void checkValidState(void) const;
  } ;


  typedef NetworkEdge<int> NetworkEdgeType;

  template<class NodeType>
    std::ostream& operator<<(std::ostream& o, const NetworkEdge<NodeType>& ne);







  class Network;
  ////////////////////////////////////////////////////////////////////
  // The NetworkPath class represents a sequence of contiguous      //
  // vertices (a path in the network). Vertices can be appended     //
  // to a path. A flow can be added to a path, and capacities of    //
  // paths can be computed.                                         //
  ////////////////////////////////////////////////////////////////////
  template<class NodeType>
    class NetworkPath {
    public:
    typedef typename std::deque<NodeType>::const_iterator vertex_iterator;

    NetworkPath(Network* network) { _network = network; }

    // Checks if the path contains no edges. 
    bool empty() const {  return ((_path.size() == 0) || (_path.size() == 1)); }

    // Removes all vertices from the path.
    void clear() { _path.clear(); }

    vertex_iterator begin() const { return _path.begin(); }

    vertex_iterator end() const { return _path.end(); }

    // Computes the maximum (residual) capacity of the path.
    double capacity(void) const;

    // appends a vertex to the end of the current flow.
    void appendVertex(NodeType nextNode);

    // prepends a vertex to the beginning of the current flow.
    void prependVertex(NodeType prevNode);

    // adds incFlow to all edges in the direction of the flow.
    void addFlow(double incFlow);


    protected:
    std::deque<NodeType> _path;
    Network*             _network;

    
    void checkValidState(void) const;
  };
  typedef NetworkPath<int> NetworkPathType;

  template<class NodeType>
    std::ostream& operator<<(std::ostream& o, const NetworkPath<NodeType>& np);











  ////////////////////////////////////////////////////////////////////
  // The Network<NumNodes> class represents a network containing    //
  // NumNodes vertices (has to be declared up front). Edges can be  //
  // added to the network using the addEdge(...) method.            //
  ////////////////////////////////////////////////////////////////////
  class Network {
    public :

    typedef std::vector<NetworkEdgeType*> EdgeList;
    typedef std::set<NetworkEdgeType*> EdgeSet;
    typedef std::set<int> NodeSet;

    // Augmenting flows can be computed by using either DFS or BFS trees.
    enum FLOW_AUGMENTING_ALG { DFS, BFS };

    // Can read DIMACS format files for testing.
    enum FILE_FORMATS { DIMACS }; 

    Network() { }

    Network(int numNodes) { init(numNodes) ; }

    Network(std::string fileName, int fileFormat = DIMACS) {  init(fileName, fileFormat);  }

    // Set the source vertex of the flow.
    void setSource(int node) {  _sourceNode = node; }

    // Set the terminal vertex of the flow.
    void setTerminal(int node) { _terminalNode = node; }

    // Returns the edge between nodeA and nodeB if one exists.
    NetworkEdgeType* edge(int nodeA, int nodeB);

    // Computes the capacity from nodeA to nodeB left in the residual network.
    double capacity(int nodeA, int nodeB);

    // Determines if there is any capacity from nodeA to nodeB in the residual network.
    bool capacityExists(int nodeA, int nodeB);

    void addFlow(int fromNode, int toNode, double incFlow);

    // Adds a directed edge to the network. 
    NetworkEdgeType* addDirectedEdge(int nodeA, int nodeB, double capacity);

    double findMaxFlow(int augAlg=BFS);

    void clearFlow(void);
    
    // Computes the set of nodes reachable from the source in the residual network.
    NodeSet& reachableSet(void) { return _reachableSet; }

    // Computes a min cut. Assumes that findMaxFlow() has already been called.
    EdgeSet& minCut(void) { return _minCut; }

  protected:
    int                                             _numNodes;
    int                                             _sourceNode;
    int                                             _terminalNode;
    int                                             _currTime;
    std::vector<int>                                _lastVisitTime;
    EdgeList                                        _edges;
    std::vector<EdgeList*>                          _outEdges;
    std::vector<EdgeList*>                          _inEdges;
    std::map<std::pair<int, int>, NetworkEdgeType*> _edgeTranslation;
    double                                          _flow;
    double                                          _cut;
    NodeSet                                         _reachableSet;
    EdgeSet                                         _minCut;
    std::vector<NetworkPathType>                    _paths;
    std::vector<int>                                _parents;

    void init(int numNodes) ;

    void init(std::string fileName, int fileFormat);

    bool findAugmentingPath(int node, NetworkPathType& augmentingPath, int augAlg);

    bool augmentingPathDFS(int node, NetworkPathType& augmentingPath);

    bool augmentingPathBFS(int node, NetworkPathType& augmentingPath);

    NodeSet& findReachableSet(int node);

    double findMinCutValue(void) ;

    bool isValidCut(const EdgeSet& e) const;

    int debugState(void) ;

    void checkValidState(void) const;

  };  
} 
#endif // GMTK_NETWORKFLOW_H

