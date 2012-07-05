/*
 * GMTK_NetworkFlow.cc
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

#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <numeric>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "GMTK_NetworkFlow.h"
#include "debug.h"
namespace networkFlow {
  const bool DebugMode = true;
  ///////////////////////////////////////////////////////////////////////////
  //                NetworkFlowException Methods                           //
  ///////////////////////////////////////////////////////////////////////////
  std::string NetworkFlowException::errorMessage(bool internal, int errNum) {
    std::string errorStr = "";
    switch(errNum) {
    case FLOW_EXCEEDS_CAPACITY :
      errorStr += "Flow should never exceed capacity : ";
      break;
    case NEGATIVE_FLOW :
      errorStr += "Negative flows are not allowed : ";
      break;
    case UNBOUNDED_FLOW :
      errorStr += "Maximum flow is unbounded : ";
      break;
    case EMPTY_PATH :
      errorStr += "Capacity of empty path not defined : ";
      break;
    case EDGE_ALREADY_EXISTS :
      errorStr += "Edge already exists. Cannot add another edge"
                  " between the same pair of vertices : ";
      break;
    case BIDIRECTIONAL_FLOW :
      errorStr += "Should not have a positive flow on both the"
                  " edge and the edge in reverse direction : ";
      break;
    case INVALID_PATH :
      errorStr += "Terminal node of edge should be the same"
	" as the start node of the next edge in"
	" a valid path : ";
      break;
    case UNSUPPORTED_FORMAT :
      errorStr += "Unsupported file format : ";
      break;
    case FILE_NOT_FOUND :
      errorStr += "Could not find File";
    default :
      errorStr += "Unknown Error : ";
    }
    if (internal) 
      errorStr += "Internal Error : ";
    return errorStr;
  }


  ///////////////////////////////////////////////////////////////////////////
  //                         NetworkEdge Methods                           //
  ///////////////////////////////////////////////////////////////////////////
  template<class NodeType> 
  double NetworkEdge<NodeType>::residualCapacity(void) const {
    if (DebugMode)
      checkValidState();
    
    if (this->hasUnlimitedCapacity())
      return capacity();
    else
      return capacity() - flow();
  }



  template<class NodeType>
  bool NetworkEdge<NodeType>::capacityExists(void) const {
    return ((hasUnlimitedCapacity()) || (residualCapacity() > 0));
  }



  template<class NodeType>
  void NetworkEdge<NodeType>::addFlow(double incFlow) {
    if (DebugMode)
      checkValidState();
    
    _flow = _flow + incFlow;
    
    if (DebugMode)
      checkValidState();
  }


  template<class NodeType>
  void NetworkEdge<NodeType>::clearFlow(void) {
    _flow = 0;
  }

  template<class NodeType>
  void NetworkEdge<NodeType>::checkValidState(void) const  {
    // Invariants :
    // (_flow >= 0)
    // (_capacity<0) || (_flow <= _capacity)
    if (_flow < 0)
      throw new NetworkFlowException(true, NetworkFlowException::NEGATIVE_FLOW);
    
    if ((_capacity >= 0) &&  (_flow > _capacity))
      throw new NetworkFlowException(true, NetworkFlowException::FLOW_EXCEEDS_CAPACITY);
  }



  template<class NodeType>
  std::ostream& operator<<(std::ostream& o, const NetworkEdge<NodeType>& ne) {
    o << "(" << ne.sourceNode() << "," << ne.targetNode() 
      << "," << ne.flow() << ":" << ne.capacity() << ")" ;
    return o;
  }





  ///////////////////////////////////////////////////////////////////////////
  //                         NetworkPath       Methods                     //
  ///////////////////////////////////////////////////////////////////////////
  template<class NodeType>
  double NetworkPath<NodeType>::capacity(void) const  {
    double pathCapacity = -1;    
    vertex_iterator  curr = begin();
    vertex_iterator  last = end();
    vertex_iterator  next = curr+1;
    
    while(next != last) {
      double edgeCapacity = _network->capacity(*curr, *next);
      if (edgeCapacity < 0) {
	// Edge has infinite capacity. Path capacity stays the same.
      } else if (pathCapacity < 0) {
	// Path so far has infinite capacity.
	pathCapacity = edgeCapacity;
      } else {
	// both path and edge have finite capacities.
	if (edgeCapacity < pathCapacity)
	  pathCapacity = edgeCapacity;
      }
      
      curr++; next++;
    }    
    return pathCapacity;
  }



  template<class NodeType>
  void NetworkPath<NodeType>::appendVertex(NodeType nextNode) {      
    _path.push_back(nextNode);
    
    if (DebugMode)
      checkValidState();
  }



  template<class NodeType>
  void NetworkPath<NodeType>::prependVertex(NodeType prevNode) {      
    _path.push_front(prevNode);
    
    if (DebugMode)
      checkValidState();
  }



  template<class NodeType>
  void NetworkPath<NodeType>::addFlow(double incFlow) {      
    vertex_iterator  curr = begin();
    vertex_iterator  last = end();
    vertex_iterator  next = curr+1;

    // Add incFlow to each edge of the network in the path.
    while(next != last) {
      _network->addFlow(*curr, *next, incFlow);
      curr++; next++;
    }
  }


  
  template<class NodeType>
  void NetworkPath<NodeType>::checkValidState(void) const {
    vertex_iterator curr = begin();
    vertex_iterator last = end();
    vertex_iterator next = curr+1;
    
    while(next != last) {
      if (!_network->edge(*curr, *next) && !_network->edge(*next, *curr)) {
	throw new NetworkFlowException(false, NetworkFlowException::INVALID_PATH);
      }
      next++; curr++;
    }
  }



  template<class NodeType>
    std::ostream& operator<<(std::ostream& o, const NetworkPath<NodeType>& np) {
    o << "start ";
    for(typename NetworkPath<NodeType>::vertex_iterator i=np.begin(); i!= np.end(); ++i) {
      const NodeType n = *i;
      o << n << " -> ";
    }
    o << "end.";
    return o;
  }








  ///////////////////////////////////////////////////////////////////////////
  //                         Network           Methods                     //
  ///////////////////////////////////////////////////////////////////////////
  void Network::init(int numNodes) {
    _numNodes = numNodes;
    _inEdges.reserve(_numNodes);
    _outEdges.reserve(_numNodes);
    for(int i=0;i<_numNodes;++i) {
      _inEdges.push_back(new std::vector<NetworkEdgeType*>);
      _outEdges.push_back(new std::vector<NetworkEdgeType*>);
      _parents.push_back(-1);
    }
    
    
    _currTime = 0;
    for(int i=0;i<_numNodes;++i) {
      _lastVisitTime.push_back(-1);
    }
  }


  // Returns a pointer to the edge from nodeA -> nodeB if one exists.
  NetworkEdgeType* Network::edge(int nodeA, int nodeB) {
    std::pair<int, int> x(nodeA, nodeB);
    if (_edgeTranslation.count(x) == 1)
      return _edgeTranslation[x];
    else
      return 0;
  }

  
  void Network::init(std::string fileName, int fileFormat) {
    if (fileFormat != DIMACS)
      throw new NetworkFlowException(false, NetworkFlowException::UNSUPPORTED_FORMAT);
    
    // Reading in the DIMACS file format. See
    // http://www.math.uni-augsburg.de/~schmidtb/bschmidt/OR_Testdata/ORTestdata.html
    // for a description of the format and for test sets.
    std::string inputLine; std::ifstream input;
    input.open(fileName.c_str(), std::ios::in);
    if (!input)
      throw new NetworkFlowException(false, NetworkFlowException::FILE_NOT_FOUND);
    
    char c; std::string problem; 
    int numNodes; int numEdges; int nodeA; int nodeB;
    double capacity;
    while(!input.eof()) {
      input >> c;
      if (input.eof())
	return;

      switch(c) {
      case 'c' : // Comment
	break;
      case 'p' : // Problem instance
	input >> problem;
	if (problem != "max")
	  throw new NetworkFlowException("Unsupported problem type : " + problem, -1);
	input >> numNodes;
	input >> numEdges;
	init(numNodes+1); 
	break;
      case 'n' : // Node descriptor.
	input >> nodeA;
	input >> c;
	if (c == 's')
	  {  setSource(nodeA); }
	else if (c == 't')
	  { setTerminal(nodeA); }
	else
	  throw new NetworkFlowException(false, NetworkFlowException::UNSUPPORTED_FORMAT);	
	break;
      case 'a' : // Edge(arc) descriptor.
	input >> nodeA;
	input >> nodeB;
	input >> capacity;
	addDirectedEdge(nodeA, nodeB, capacity);
	break;
      default :  // Unknown command.
	throw new NetworkFlowException(false, NetworkFlowException::UNSUPPORTED_FORMAT);
      }
    }
    input.close();
  }



  // Computes the (direct) capacity from nodeA to nodeB in the residual network.
  double Network::capacity(int nodeA, int nodeB) {
    NetworkEdgeType* fwdEdge = edge(nodeA, nodeB);
    NetworkEdgeType* revEdge = edge(nodeB, nodeA);
    double capacity = 0.0;
    
    // If there is an edge from nodeA to nodeB, we can use
    // the residual capacity of this edge.
    if (fwdEdge) {
      if (fwdEdge->hasUnlimitedCapacity())
	return -1;
      else
	capacity += fwdEdge->residualCapacity();
    }
    
    // If there is an edge from nodeB to nodeA, we can use
    // the flow on this edge.
    if (revEdge) {
      capacity += revEdge->flow();
    }
      
    return capacity;
  }



  bool Network::capacityExists(int nodeA, int nodeB) {
    double cap = capacity(nodeA, nodeB);
    // capacity < 0 is used to indicate infinite capacity.
    return ((cap < 0) || (cap > 0));
  }


  
  void Network::addFlow(int fromNode, int toNode, double incFlow) {
    if (incFlow < 0)
      throw new NetworkFlowException(true, NetworkFlowException::NEGATIVE_FLOW);
    
    NetworkEdgeType* fwdEdge = edge(fromNode, toNode);
    NetworkEdgeType* revEdge = edge(toNode, fromNode);

    // First try and stick the flow on the reverse edge if one exists.
    if (revEdge) {
      double revFlow = revEdge->flow();
      if (revFlow >= incFlow) {
	revEdge->addFlow(-incFlow);
	return;
      } else if (revFlow > 0) {
	incFlow -= revFlow;
	revEdge->addFlow(-revFlow);
      }
    }


    // If there is a forward edge, stick all the remaining flow on the forward edge. 
    if (fwdEdge) {
      if (fwdEdge->hasUnlimitedCapacity()) {
	fwdEdge->addFlow(incFlow);
	return;
      } else {
	double fwdCap = fwdEdge->residualCapacity();
	if (fwdCap >= incFlow) {
	  fwdEdge->addFlow(incFlow);
	  return;
	} else {
	  throw new NetworkFlowException(true, NetworkFlowException::FLOW_EXCEEDS_CAPACITY);
	}
      }
    }
  }




  // Adds a directed edge from nodeA to nodeB if there is no edge from
  // nodeA to nodeB already in the network. 

  // TODO: Instead of throwing an exception in the case when an edge already 
  // exists should I just add to the capacity of the current edge?
  NetworkEdgeType* Network::addDirectedEdge(int nodeA, int nodeB, double capacity) {

    if (this->edge(nodeA, nodeB) != 0) {
      throw new NetworkFlowException(false, NetworkFlowException::EDGE_ALREADY_EXISTS);
    }

    
    NetworkEdgeType* e = new NetworkEdgeType(nodeA, nodeB, capacity);
     std::pair<int, int> x(nodeA, nodeB);
     _edgeTranslation[x] = e;
    _edges.push_back(e);
    _outEdges[nodeA]->push_back(e);
    _inEdges[nodeB]->push_back(e);

    return e;
  }


  void Network::clearFlow(void) {
    throw new NetworkFlowException(" Not yet implemented : ", -1);
  }


  // findMaxFlow(int augAlg) computes the maximum flow from _sourceNode to
  // _terminalNode in the network using the Ford-Fulkerson algorithm.
  // A good description of this algorithm can be found in 
  // Cormen,Leiserson,Rivest (Introduction to algorithms)
  // Flow augmentations can be computed by using either a DFS search or
  // a BFS search. The BFS search is the "right" way of doing this, but
  // I've added in the DFS search as well for debugging reasons.
  double Network::findMaxFlow(int augAlg) {
    // _currTime and _lastVisitTime are used to determine which nodes
    // have been visited in the current iteration and which have not.
    // In general, _lastVisitTime[i] is the last iteration at which 
    // node i was visited. So, if _lastVisitTime[i] < _currTime, then
    // node i was not visited in this iteration. 
    _currTime = 0;
    NetworkPathType path(this);
    for(int i=0;i<_numNodes;++i)
      _lastVisitTime[i] = -1;
    _flow = 0;

    bool augmentingPathExists = findAugmentingPath(_sourceNode, path, augAlg);
    // A flow is maximal if and only if there is no augmenting path in the
    // residual network. So check if an augmenting path exists. 
    while(augmentingPathExists) {
      infoMsg(IM::Mod, "Found augmenting path\n");
      
      // Stick the maximum possible flow along this augmenting path. 
      double incFlow = path.capacity();
      path.addFlow(incFlow);       
      _flow += incFlow;

      // Store the set of all augmenting paths so that other routines 
      // can use them (for example to compute internally vertex disjoint
      // paths).
      _paths.push_back(path);
      path.clear();
      _currTime++;
      augmentingPathExists = findAugmentingPath(_sourceNode, path, augAlg);
    }

    ++_currTime;
    infoMsg(IM::Mod, "Computing reachable and cut sets\n");


    _cut = findMinCutValue();
    if (_cut != _flow) {
      infoMsg(IM::High, "Error: Flow does not match cut. ");
      EdgeSet::iterator ei,eend;
      ei=_minCut.begin(); eend=_minCut.end();
      for(;ei!=eend;++ei) {
	infoMsg(IM::Mod,  "Edge %d -> %d", (*ei)->sourceNode(), (*ei)->targetNode());
      }
      // debugState();
      throw new NetworkFlowException("Internal Error: Flow does not match cut", -1);
    } else {
      return _flow;
    }
  }



  bool Network::findAugmentingPath(int node, NetworkPathType& augmentingPath, int augAlg) {
    if (augAlg == DFS)
      return augmentingPathDFS(node, augmentingPath);
    else
      return augmentingPathBFS(node, augmentingPath);
  }



  // Compute an augmenting path along a DFS tree. 
  // This code is a simple recursive implementation of DFS (see LCR, Intro to algorithms).
  bool Network::augmentingPathDFS(int node, NetworkPathType& augmentingPath) {
    if (_lastVisitTime[node] == _currTime) {
      return false;
    }else if (node == _terminalNode) {
      augmentingPath.prependVertex(node);
      return true;
    } else {
      _lastVisitTime[node] = _currTime;

      // Need to check outgoing edges which are not saturated and 
      // incoming edges which have some flow on them (because these
      // are the edges in the residual network). 
      EdgeList* outEdges = _outEdges[node];
      EdgeList::iterator i,begin,end;
      begin = outEdges->begin(); end = outEdges->end();
      for(i=begin;i!=end;++i) {
	int otherNode = (*i)->targetNode();
	if ((_lastVisitTime[otherNode] < _currTime) &&
	    (capacityExists(node, otherNode)) &&
	    (augmentingPathDFS(otherNode, augmentingPath))) {
	  augmentingPath.prependVertex(node);
	  return true;
	}
      }


      EdgeList* inEdges = _inEdges[node];
      begin = inEdges->begin(); end = inEdges->end();
      for(i=begin;i!=end;++i) {
	int otherNode = (*i)->sourceNode();
	if ((_lastVisitTime[otherNode] < _currTime) &&
	    (capacityExists(node, otherNode)) &&
	    augmentingPathDFS(otherNode, augmentingPath)) {
	  augmentingPath.prependVertex(node);
	  return true; 
	}
      }

      return false;
    }

  }


  // Compute an augmenting path along a DFS tree. 
  // This code is a simple implementation of DFS (see LCR, Intro to algorithms).
  bool Network::augmentingPathBFS(int node, NetworkPathType& augmentingPath) {
    std::queue<int> unexploredNodes;
    std::map<int, int> _parent;
    unexploredNodes.push(node);
    _parent[node] = node;
    _lastVisitTime[node] = _currTime;
    
    while(unexploredNodes.size() > 0) {
      int currNode = unexploredNodes.front();
      if (currNode == _terminalNode) {
	while(currNode != _parent[currNode]) {
	  augmentingPath.prependVertex(currNode);
	  currNode = _parent[currNode];
	}
	augmentingPath.prependVertex(currNode);
	return true;
      } else {
	EdgeList* outEdges = _outEdges[currNode];
	EdgeList::iterator i,begin,end;
	begin = outEdges->begin(); end = outEdges->end();
	for(i=begin;i!=end;++i) {
	  int nextNode = (*i)->targetNode();
	  if ((_lastVisitTime[nextNode] < _currTime) &&
	      (capacityExists(currNode, nextNode)) && 
	      (_parent.count(nextNode) == 0)) {
	    _parent[nextNode] = currNode;
	    _lastVisitTime[nextNode] = _currTime;
	    unexploredNodes.push(nextNode);
	  }
	}


	EdgeList* inEdges = _inEdges[currNode];
	begin = inEdges->begin(); end = inEdges->end();
	for(i=begin;i!=end;++i) {
	  int nextNode = (*i)->sourceNode();
	  if ((_lastVisitTime[nextNode] < _currTime) &&
	      (capacityExists(currNode, nextNode)) && 
	      (_parent.count(nextNode) == 0)) {
	    _lastVisitTime[nextNode] = _currTime;
	    _parent[nextNode] = currNode;
	    unexploredNodes.push(nextNode);
	  }
	}

	_lastVisitTime[node] = _currTime;
      }
      unexploredNodes.pop();
    }


    return false;
  }




  // The reachable set is the set of vertices v so that a path exists
  // from _sourceNode to v in the residual network. 
  std::set<int>& Network::findReachableSet(int node) {
    _lastVisitTime[node] = _currTime;
    _reachableSet.insert(node);

    
    EdgeList::iterator i,begin,end;
    EdgeList* edges = _outEdges[node];
    begin = edges->begin(); end = edges->end();    
    for(i=begin;i!=end;++i) {
      int nextNode = (*i)->targetNode();
      if ((_lastVisitTime[nextNode] < _currTime) &&
	  capacityExists(node, nextNode)) {
	_parents[nextNode] = node;
	findReachableSet(nextNode);
      }
    }


    edges = _inEdges[node];
    begin = edges->begin(); end = edges->end();    
    for(i=begin;i!=end;++i) {
      int nextNode = (*i)->sourceNode();
      if ((_lastVisitTime[nextNode] < _currTime) &&
	  capacityExists(node, nextNode)) {
	_parents[nextNode] = node;
	findReachableSet(nextNode);
      }
    }

    
    return _reachableSet;
  }

  

  // Once a max-flow is constructed, there is no augmenting path
  // in the network. Then the cut defined by $(S,T)$ where
  // $S$ is the set of vertices reachable from _sourceNode in the
  // residual network and $T$ is the complement of $S$ forms a min cut.
  double Network::findMinCutValue(void) {
    _minCut.clear(); _reachableSet.clear(); 
    _parents[0] = 0;
    _reachableSet.clear();
    findReachableSet(_sourceNode);
    std::set<int>::iterator ni,nbegin,nend;
    _cut = 0;
    nbegin=_reachableSet.begin(); nend=_reachableSet.end();

    for(ni=nbegin;ni!=nend;++ni) {
      int currNode = *ni;
      EdgeList* edges = _outEdges[currNode];
      EdgeList::iterator ei=edges->begin(); 
      EdgeList::iterator eend=edges->end();
      for(;ei!=eend;++ei) {
	int nextNode = (*ei)->targetNode(); 
	if ((_reachableSet.find(nextNode) == nend)) {
	  double cap = (*ei)->capacity();
	  if (cap >= 0) {
	    _cut += cap;
	    _minCut.insert(*ei);
	  } else {
	    throw new NetworkFlowException("Internal Error: Infinite capacity edge cut", -1);
	  }
	}
      }
    }
    return _cut;
  }


  bool Network::isValidCut(const Network::EdgeSet&  e) const {
    // try and find a path from the source to the destination - not yet implemented.
    
    return false;
  }


  // This is routine that lets me see debug the network without going through
  // gdb. 
  int Network::debugState(void) {
    std::string command; int arg0, arg1;
    while (command != "quit") {
      std::cin >> command;
      if (command == "oedges") {
	// Display all edges going out of arg0.
	std::cin >> arg0;
	std::cout << "Edges leaving " << arg0 << std::endl;
	EdgeList* edges = _outEdges[arg0];
	EdgeList::iterator ei=edges->begin(); 
	EdgeList::iterator eend=edges->end();
	for(;ei!=eend;++ei) {
	  int nextNode = (*ei)->targetNode(); 
	  std::cout << "  : "  << nextNode << std::endl;
	}
      } else if (command == "iedges") {
	// Display all edges going in to arg0.
	std::cin >> arg0;
	std::cout << "Edges going to " << arg0 << std::endl;
	EdgeList* edges = _inEdges[arg0];
	EdgeList::iterator ei=edges->begin(); 
	EdgeList::iterator eend=edges->end();
	for(;ei!=eend;++ei) {
	  int nextNode = (*ei)->sourceNode(); 
	  std::cout << "  : "  << nextNode << std::endl;
	}
      } else if (command == "pathto") {
	// Display a path to arg0 in residual network
	std::cin >> arg0;
      } else if (command == "pathfrom") {
	// Display a path from arg0 to terminal in residual network.
	std::cin >> arg1;
      } else if (command == "capacity") {
	// Display the capacity in the residual network from arg0 to arg1.
	std::cin >> arg0;
	std::cin >> arg1;
	std::cout << "Capacity from " << arg0 << " to " << arg1 << " = " << capacity(arg0, arg1) << std::endl;
      }
    }
    return 1;

  }


  void Network::checkValidState(void) const {
    // Add checking code here.
  }
} 






