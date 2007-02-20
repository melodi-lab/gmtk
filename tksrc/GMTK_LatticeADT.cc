/*-
 * GMTK_LatticeADT.cc
 *
 *  Written by Gang Ji <gang@ee.washington.edu>
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


#include <cstring>
#include <cstdlib>
#include <cmath>

#include "GMTK_LatticeADT.h"
#include "GMTK_NamedObject.h"
#include "GMTK_Vocab.h"
#include "GMTK_GMParms.h"
#include "error.h"


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::LatticeADT(ifs, vocab)
 *     default constructor
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
LatticeADT::LatticeADT()
  : _latticeNodes(NULL), _numberOfNodes(0),
    _numberOfLinks(0), _start(0), _end(0),
    _lmscale(0), _wdpenalty(0), _acscale(0),
    _amscale(0), _base(0), _frameRelax(0),
    _latticeFile(NULL), _numLattices(0), _curNum(0),
    _nodeCardinality(0), _wordCardinality(0), _timeCardinality(0) {
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::~LatticeADT(ifs, vocab)
 *     default destructor
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
LatticeADT::~LatticeADT() {
  if  ( _latticeNodes )
    delete [] _latticeNodes;
  if ( _latticeFile )
    delete _latticeFile;
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::readFromHTKLattice(ifs, vocab)
 *     read an HTK SLF lattice file format file. This routine *should* conform
 *     to the description that is listed in Chapter 20 (in 3.2.1) of the HTKbook.
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
void LatticeADT::readFromHTKLattice(iDataStreamFile &ifs, const Vocab &vocab) 
{

  /////////////////////////
  // Read Lattice Header
  /////////////////////////

  char *s_tmp, *ptr;
  do {
    ifs.read(s_tmp);
    if ( s_tmp == NULL )
      break;

    // TODO: change below atoi's to strtouls, ideally add a common routine in general.{h,cc}

    // TODO: go through chapter 20 in HTK book nd finish implementing
    // this function to fully support HTK lattices.

    // TODO: allow case to be both in checks below. 

    if ( strstr(s_tmp, "VERSION") != 0 ) {
      // skip the line
    } else if ( strstr(s_tmp, "UTTERANCE") != NULL ) {
      // skip the line
    } else if ( strstr(s_tmp, "lmscale=") != NULL ) {
      // parsing language model scale
      ptr = strchr(s_tmp, '=');
      _lmscale = atof(++ptr);
    } else if ( strstr(s_tmp, "wdpenalty=") != NULL ) {
      // parsing word penalty
      ptr = strchr(s_tmp, '=');
      _wdpenalty = atof(++ptr);
    } else if ( strstr(s_tmp, "acscale") != NULL ) {
      // parse acoustic scale
      ptr = strchr(s_tmp, '=');
      _acscale = atof(++ptr);
    } else if ( strstr(s_tmp, "base=") != NULL ) {
      // parse log base
      ptr = strchr(s_tmp, '=');
      _base = atof(++ptr);
    } else if ( strstr(s_tmp, "acscale") != NULL ) {
      // parse acoustic likelyhood scale
      ptr = strchr(s_tmp, '=');
      _amscale = atof(++ptr);
    } else if ( strstr(s_tmp, "START=") != NULL || strstr(s_tmp, "start=") != NULL ) {
      // parse start node id
      ptr = strchr(s_tmp, '=');
      _start = (unsigned)atoi(++ptr);
    } else if ( strstr(s_tmp, "END=") != NULL || strstr(s_tmp, "end=") != NULL ) {
      // parse start node id
      ptr = strchr(s_tmp, '=');
      _end = (unsigned)atoi(++ptr);
    } else if ( (strstr(s_tmp, "N=") != NULL ) || (strstr(s_tmp, "NODES=") != NULL )) {
      // parse number of nodes
      ptr = strchr(s_tmp, '=');
      _numberOfNodes = (unsigned)atoi(++ptr);
      if ( _numberOfNodes > _nodeCardinality )
	error("Error in lattice '%s', number of nodes in lattice %d is bigger than RV lattice node cardinality %d",
	      ifs.fileName(),_numberOfNodes, _nodeCardinality);
    } else if ( (strstr(s_tmp, "L=") != NULL) || (strstr(s_tmp, "LINKS=") != NULL) ) { 
      // parse number of links/edges
      ptr = strchr(s_tmp, '=');
      _numberOfLinks = (unsigned)atoi(++ptr);
      delete [] s_tmp;
      break;		// read to parse node and link information
    }
		
    delete [] s_tmp;
    s_tmp = NULL;
  } while ( ! ifs.isEOF() );

  if (_numberOfNodes < 2) 
    error("Error in lattice '%s', lattice only has %d nodes, must have at least 2 nodes",
	  ifs.fileName(),_numberOfNodes);

  if (_numberOfLinks < 1) 
    error("Error in lattice '%s', lattice only has %d edges, must have at least 1 edge",
	  ifs.fileName(),_numberOfLinks);

  if (_start == _end) {
    error("Error in lattice '%s', start node id = %d is same as end node id = %d",
	  ifs.fileName(),_start,_end);
  }


  /////////////////////////
  // Read Lattice Nodes
  /////////////////////////

  if ( _latticeNodes )
    delete [] _latticeNodes;
  _latticeNodes = new LatticeNode [_numberOfNodes];

  const unsigned linesize = 4096;
  char line[linesize];
  const char seps[] = " \t\n";
  unsigned id;

  // to assert there is no errors in the latice, we perform two checks:
  // 1. make sure that for each edge, the end time is later than start time
  // 2. make sure that end node in the lattice has the very latest time
  float latestTime = 0.0;

  // reading nodes
  for ( unsigned i = 0; i < _numberOfNodes; i++ ) {

    ifs.readLine( line, linesize);
    ptr = strtok(line, seps);

    if ( ptr[0] != 'I' || ptr[1] != '=' )
      error("Error in lattice '%s', line %d, expecting I= within line, but got '%s'", 
	    ifs.fileName(),ifs.lineNo(),ptr);
    if ( (id = atoi(ptr+2)) > _numberOfNodes )
      error("Error in lattice '%s', line %d, node id %d is bigger than number of nodes %d",
	    ifs.fileName(),ifs.lineNo(),
	    id, _numberOfNodes);
    if ( (ptr = strtok(NULL, seps)) != NULL ) {
      // there is time information
      if ( ptr[0] == 't' ) {
	// set up the time
	_latticeNodes[id].time = atof(ptr+2);

	// by default, there is no frame constrain
	_latticeNodes[id].startFrame = 0;
	_latticeNodes[id].endFrame = ~0;

	// update latest time
	if ( latestTime < _latticeNodes[id].time )
	  latestTime = _latticeNodes[id].time;
      }
    }
  }

  // make sure end node has the latest time
  if ( _latticeNodes[_end].time < latestTime )
    error("When reading lattice '%s', lattice end node %d has a time of %f which is not the latest found time which is %f.",ifs.fileName(),_end,_latticeNodes[_end].time,latestTime);


  ////////////////////////////
  // Read Lattice Edges/Links
  ////////////////////////////

  // reading links
  double score;
  unsigned endNodeId = 0;
  for ( unsigned i = 0; i < _numberOfLinks; i++ ) {
    ifs.readLine(line, linesize);
    ptr = strtok(line, seps);
    if ( ptr[0] != 'J' || ptr[1] != '=' )
      error("Error in lattice %s, line %d, expecting J= within line '%s'",
	    ifs.fileName(),ifs.lineNo(),ptr);
    if ( (id = atoi(ptr+2)) > _numberOfLinks )
      error("Error in lattice %s, line %d, link id %d is bigger than number of links %d", 
	    ifs.fileName(),ifs.lineNo(),id, _numberOfLinks);

    // starting node
    if ( (ptr = strtok(NULL, seps)) == NULL )
      error("Error in lattice '%s', line %d, expecting starting node for link %d",
	    ifs.fileName(),ifs.lineNo(),id);
    if ( ptr[0] != 'S' || ptr[1] != '=' )
      error("Error in lattice '%s', line %d, expecting starting node for link %d", 
	    ifs.fileName(),ifs.lineNo(),id);
    if ( (id = atoi(ptr+2)) > _numberOfNodes )
      error("Error in lattice '%s', line %d,  node id %d is bigger than number of lattice nodes %d", 
	    ifs.fileName(), ifs.lineNo(), id, _numberOfNodes);

    LatticeEdge edge;
    while ( (ptr = strtok(NULL, seps)) != NULL ) {
      switch ( ptr[0] ) {
      case 'E':
	// ending link
	endNodeId = atoi(ptr+2);
	break;
      case 'W':
	// word token
	edge.emissionId = vocab.index(ptr+2);
	if ( edge.emissionId == vocab.index("<unk>") ) {
	  // ptr is the string value of an edge in the lattice. This
	  // case corresponds to the case where the lattice contains
	  // an edge that is labeled '<unk>' or it contains an edge
	  // that has a label that is unknown in the current
	  // vocabulary. TODO: ultimately, give options to have unk,
	  // no unk, etc. similar to SRILM, but for now map either the
	  // string <unk> or any other unknown word relative to the
	  // vocab object to the id corresponding to the <unk> string.
	  // Note: an LM, if it is used for this, must have an <unk> ability
	  //  (i.e., an LM must be able to re-score an unk symbol).
	  // Note: this will lose the original ID of the word (i.e., any
	  // printing will print out unk).
	  // TODO: fix this /rethink.
	  // error("Error: word '%s' in lattice '%s' cannot be found in vocab", ptr+2, ifs.fileName());
	}
	break;
      case 'a':
	// acoustic score
	score = atof(ptr+2);
	if ( _base == 0 ) {
	  score = log(score);
	} else {
	  score *= log(_base);
	}
	edge.ac_score.setFromLogP(score);
	break;
      case 'l':
	// lm score
	score = atof(ptr+2);
	if ( _base == 0 ) {
	  score = log(score);
	} else {
	  score *= log(_base);
	}
	edge.lm_score.setFromLogP(score);
	break;
      case 'd':
	// duration
	break;
      case 'p':
	// posterior
	score = atof(ptr+2);
	edge.posterior.setFromP(score);
	break;
      default: 
	// is its anything other than the above, simply ignore it for now. 
	infoMsg(IM::Warning,"WARNING: Lattice '%s', line %d, contains a token '%s' that is not currently supported and will be ignored\n",ifs.fileName(),ifs.lineNo(),ptr);
	//	error("Error in lattice '%s', line %d, found unknown token '%s'",ifs.fileName(),ifs.lineNo(),ptr);
	break;
      }
    }


    infoMsg(IM::Giga, "Edge information: emissionID=%d, ac_score=%f, lm_score=%f, posterior=%f, node_id=%d, end_node_id=%d\n",edge.emissionId,edge.ac_score.val(),edge.lm_score.val(),edge.posterior.val(),id,endNodeId);


    // make sure the end time is later than start time
    if ( _latticeNodes[id].time >= _latticeNodes[endNodeId].time )
      error("Error in lattice '%s', lattice edge %d, start node %d time (%f) should be earlier than end node %d time(%f)",
	    ifs.fileName(),
	    i, 
	    id, _latticeNodes[id].time, 
	    endNodeId, _latticeNodes[endNodeId].time);

    _latticeNodes[id].edges.insert(endNodeId, edge);

  }

  if ( _latticeNodes[_start].time >= _latticeNodes[_end].time )
    error("Error in lattice '%s', overall start node %d time (%f) should be earlier than end node %d time(%f)",
	  ifs.fileName(),
	  _start, _latticeNodes[_start].time, 
	  _end, _latticeNodes[_end].time);

  normalizePosterior();
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::read(ifs)
 *     read the lattice object from master file
 *
 * Results:
 *     no results:
 *
 * Note:
 *     Be careful! the frame index information is NOT set here.
 *     This will wait until the global observation matrix is loaded
 *     so that the last frame index is known.
 *-----------------------------------------------------------------------
 */
void LatticeADT::read(iDataStreamFile &is) {
  // read in name
  NamedObject::read(is);

  // read in node cardinality or lattice list filename
  is.read(_latticeFileName, "Can't read lattice file or number of fetures");

  // if lattice file name is a string, then it is a list of lattices
  if ( ! strIsInt(_latticeFileName.c_str(), (int*)&_nodeCardinality) ) {
    if ( iterable() )
      error("ERROR: in lattices named '%s' in file '%s' line %d, can't have lattices defined recursively in files", name().c_str(),is.fileName(),is.lineNo());

   initializeIterableLattice(_latticeFileName);
  } else {
    // read in next to see whether it is time parent
    string tmpStr;
    int tmpInt;
    is.read(tmpStr);
    char *latticeFile;
    if ( strIsInt(tmpStr.c_str(), &tmpInt) ) {
      _timeCardinality = (unsigned) tmpInt;
      // read in lattice filename
      is.read(latticeFile, "Can't read lattice filename");
    } else {
      // no time parent
      latticeFile = new char [tmpStr.size() + 1];
      strcpy(latticeFile, tmpStr.c_str());
    }

    // read in vocab filename
    string vocabName;
    is.read(vocabName, "Can't read vocab name");

    if ( GM_Parms.vocabsMap.find(vocabName) == GM_Parms.vocabsMap.end() )
      error("Error: reading file '%s' line '%d', LatticeCPT '%s' specifies Vocab name '%s' that does not exist", is.fileName(), is.lineNo(), _name.c_str(), vocabName.c_str());

    // set vocabulary and read in HTK lattice
    Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
    iDataStreamFile lfifs(latticeFile,  // name
			  false,        // isBinary
			  false,        // use CPP if ASCII
			  NULL,         // extra CPP options
			  '#'           // extra comment character to use, use SLF comment '#'
			  );

    readFromHTKLattice(lfifs, *vocab);
    delete [] latticeFile;

    _wordCardinality = vocab->size();

    // read in options for choosing scores
    // 0x0,0x3 no AC
    // 0x1 AC
    // 0x2 ACScale
    // 0x0,0x3 (<<2) no LM
    // 0x1 (<<2) LM
    // 0x2 (<<2) LMScale
    // 0x0 (<<4) no WDPenalty
    // 0x1 (<<4) WDPenalty
    // 0x0 (<<5) no posterior
    // 0x1 (<<5) posterior
    // note: if posterior is used, all the others are disabled
    unsigned option = 0;
    string trigger = "UseScore";
    if ( is.readIfMatch(trigger, "reading score options") ) {
      char *str = new char[1024];
      is.read(str, "reading score options");
      if ( strcmp(str, "Posterior") == 0 ) {
	option = 0x1u << 5;
      } else {
	char *tok = strtok(str, "+");
	while ( tok != NULL ) {
	  if ( strcmp(tok, "AC") == 0 ) {
	    option &= ~0x3;
	    option |= 0x1;
	  } else if ( strcmp(tok, "ACScale") == 0 ) {
	    option &= ~0x3;
	    option |= 0x2;
	  } else if ( strcmp(tok, "LM") == 0 ) {
	    option &= ~0xc;
	    option |= 0x4;
	  } else if ( strcmp(tok, "LMScale") == 0 ) {
	    option &= ~0xc;
	    option |= 0x8;
	  } else if ( strcmp(tok, "WDPenalty") == 0 ) {
	    option |= 0x10;
	  } else
	    error("Error: reading file '%s' line '%d', unknown score option '%s'", is.fileName(), is.lineNo(), tok);
	  tok = strtok(NULL, "+");
	}
      }
      delete [] str;
    }
    useScore(option);

    // read in frame relaxation
    is.read(_frameRelax, "Can't read lattice frame relaxation");
  }

  // Note: the frame index information is NOT set here.
  // This will wait until the global observation matrix is loaded
  // so that the last frame index is known.
}


/*-
 *-----------------------------------------------------------------------
 * seek
 *      Seek to a particular lattice
 *
 * Results:
 *      none
 *-----------------------------------------------------------------------
 */
void LatticeADT::seek(unsigned nmbr) 
{
  string tmpStr;
  int tmpInt;
  char trigger[] = "UseScore";

  if (!iterable())
    error("ERROR: trying to seek in non-iterable lattice, '%s'\n",
	  _curName.c_str() );

  // rewind the file
  _latticeFile->rewind();
  _latticeFile->read(_numLattices, "num lattices");
  if ( _numLattices <= nmbr )
    error("number of lattices (%d) in '%s' is less than segment number (%d)", _numLattices, _latticeFileName.c_str(), nmbr);

  // Iterate through the lattices based on  the nubmer of lines each CPT has
  for ( unsigned i = 0; i < nmbr; i++ ) {
    // Read index, lattice name, and cardinality 
    for ( unsigned j = 0; j < 3; j++ ) {
      _latticeFile->read(tmpStr);
    }

    // Read the next line, if is an integer there is a time parent and one 
    // more line needs to be read.  
    _latticeFile->read(tmpStr);
    if ( strIsInt(tmpStr.c_str(), &tmpInt) ) {
      _latticeFile->read(tmpStr);
    } 

    // Read Vocabulary 
    _latticeFile->read(tmpStr);

    // Read UseScore 
    _latticeFile->read(tmpStr);
    if (strcmp(tmpStr.c_str(), trigger) != 0) {
      error("ERROR: Error seeking to lattice %d\n", nmbr);
    }

    // Read UseScore options and slack
    for ( unsigned j = 0; j < 2; j++ ) {
      _latticeFile->read(tmpStr);
    }
  }

  // set current nuber to this so that next can be called
  _curNum = nmbr - 1;
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::initializeIterableLattice(fileName)
 *     initialize an iterable lattice cpt
 *
 * Results:
 *     None.
 *-----------------------------------------------------------------------
 */
void LatticeADT::initializeIterableLattice(const string &fileName) {
  if ( iterable() )
    error("ERROR: can't call initializeIterableDT() for recursively");

  _latticeFileName = fileName;
  _latticeFile = new iDataStreamFile(_latticeFileName.c_str(), false, false);
  _numLattices = 0;
  _curNum = -1;

  beginIterableLattice();
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::beginIterableLattice()
 *     initialize an iterable lattice cpt
 *
 * Results:
 *     None.
 *-----------------------------------------------------------------------
 */
void LatticeADT::beginIterableLattice() {
  if ( ! iterable() )
    error("ERROR: can't call beginIterableDT() for non-file lattice");
  _latticeFile->rewind();
  _latticeFile->read(_numLattices, "num lattices");
  _curNum = -1;

  // read in the fist lattice
  nextIterableLattice();
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::nextIterableLattice()
 *     read an lattice cpt
 *
 * Results:
 *     None.
 *-----------------------------------------------------------------------
 */
void LatticeADT::nextIterableLattice() 
{
  if ( ! iterable() )
    error("ERROR: can't call nextIterableDT() for non-file lattice");

  int readLatticeNum;

  // increment the index
  ++_curNum;
  _latticeFile->read(readLatticeNum, "lattice index");
  if ( _curNum != readLatticeNum )
    error("ERROR: reading from file '%s', expecting DT number %d but got number %d\n", _latticeFileName.c_str(), _curNum, readLatticeNum);

  _latticeFile->read(_curName, "current lattice name");
  _latticeFile->read(_nodeCardinality, "number of nodes");

  // read in next to see whether it is time parent
  string tmpStr;
  int tmpInt;
  char *latticeFile;
  _latticeFile->read(tmpStr);
  if ( strIsInt(tmpStr.c_str(), &tmpInt) ) {
    _timeCardinality = (unsigned) tmpInt;
    // read in lattice filename
    _latticeFile->read(latticeFile, "Can't read lattice filename");
  } else {
    // no time parent
    latticeFile = new char [tmpStr.size() + 1];
    strcpy(latticeFile, tmpStr.c_str());
  }

  // read in vocab filename
  string vocabName;
  _latticeFile->read(vocabName, "Can't read vocab name");

  if ( GM_Parms.vocabsMap.find(vocabName) == GM_Parms.vocabsMap.end() )
    error("Error: reading file '%s' line '%d', LatticeCPT '%s' specifies Vocab name '%s' that does not exist", _latticeFile->fileName(), _latticeFile->lineNo(), _name.c_str(), vocabName.c_str());

  // set vocabulary and read in HTK lattice
  Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
  iDataStreamFile lfifs(latticeFile, false, false);
  readFromHTKLattice(lfifs, *vocab);
  _wordCardinality = vocab->size();
  delete [] latticeFile;

  // read in options for choosing scores
  // 0x0,0x3 no AC
  // 0x1 AC
  // 0x2 ACScale
  // 0x0,0x3 (<<2) no LM
  // 0x1 (<<2) LM
  // 0x2 (<<2) LMScale
  // 0x0 (<<4) no WDPenalty
  // 0x1 (<<4) WDPenalty
  // 0x0 (<<5) no posterior
  // 0x1 (<<5) posterior
  // note: if posterior is used, all the others are disabled
  unsigned option = 0;
  string trigger = "UseScore";
  if ( _latticeFile->readIfMatch(trigger, "reading score options") ) {
    char* str = NULL;
    // note, read will allocate the string, but we need to free it below.
    // since this is called only once, not too bad to allocate and delete here.
    _latticeFile->read(str, "reading score options");

    if ( strcmp(str, "Posterior") == 0 ) {
      option = 0x1u << 5;
    } else {
      char *tok = strtok(str, "+");
      while ( tok != NULL ) {
	if ( strcmp(tok, "AC") == 0 ) {
	  option &= ~0x3;
	  option |= 0x1;
	} else if ( strcmp(tok, "ACScale") == 0 ) {
	  option &= ~0x3;
	  option |= 0x2;
	} else if ( strcmp(tok, "LM") == 0 ) {
	  option &= ~0xc;
	  option |= 0x4;
	} else if ( strcmp(tok, "LMScale") == 0 ) {
	  option &= ~0xc;
	  option |= 0x8;
	} else if ( strcmp(tok, "WDPenalty") == 0 ) {
	  option |= 0x10;
	} else
	  error("Error: reading file '%s' line '%d', unknown score option '%s'", _latticeFile->fileName(), _latticeFile->lineNo(), tok);
	tok = strtok(NULL, "+");
      }
    }
    free(str);
  }
  useScore(option);

  // read in frame relaxation
  _latticeFile->read(_frameRelax, "Can't read lattice frame relaxation");
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::resetFrameIndices(lastFrameId)
 *     reset all frame indices when global observation matrix is loaded
 *     into memory.  This is a simple linear warping.
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
void LatticeADT::resetFrameIndices(unsigned numFrames) {
  const unsigned lastFrameId = numFrames - 1;
  float warpRate = lastFrameId / _latticeNodes[_end].time;

  unsigned frame;
  for ( unsigned i = 0; i < _numberOfNodes; i++ ) {
    frame = (unsigned)roundf(warpRate * _latticeNodes[i].time);
    _latticeNodes[i].startFrame = (frame < _frameRelax) ? 0 : frame - _frameRelax;

    frame += _frameRelax;
    _latticeNodes[i].endFrame = (frame > lastFrameId - 1) ? lastFrameId - 1: frame;
  }

  // special treatment for start
  _latticeNodes[_start].startFrame = _latticeNodes[_start].endFrame = 0;

  // special treatment for end
  _latticeNodes[_end].startFrame = _latticeNodes[_end].endFrame = lastFrameId;
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::useScore
 *     use difference score options
 *     refer read for what option means
 *-----------------------------------------------------------------------
 */
void LatticeADT::useScore(unsigned option) 
{
  for ( unsigned i = 0; i < _numberOfNodes; ++i ) {
    if ( _latticeNodes[i].edges.totalNumberEntries() > 0 ) {
      shash_map_iter<unsigned, LatticeEdge>::iterator it;
      _latticeNodes[i].edges.begin(it);
      do {
	// score of the edge in ln value
	double score = 0;

	if ( option & (0x1u<<5) ) {
	  // use only posterior
	  score = (*it).posterior.val();
	} else {
	  // do we use AC score?
	  unsigned x = option & 0x3;
	  switch ( x ) {
	  case 1:	// AM score only
	    score = (*it).ac_score.val();
	    break;
	  case 2: // AM^a
	    score = (*it).ac_score.val() * _acscale;
	    break;
	  default: // no AM score
	    break;
	  }

	  // do we use LM score?
	  x = (option >> 2) & 0x3;
	  switch ( x ) {
	  case 1: // LM score only
	    score += (*it).lm_score.val();
	    break;
	  case 2: // LM^b
	    score += (*it).lm_score.val() * _lmscale;
	    break;
	  default: // no LM score
	    break;
	  }

	  // do we use insertion penalty?
	  if ( option & 0x10 )
	    score += _wdpenalty;
	}

	// check whether the score is too small that will have zero
	// probability.
	if ( score <= LSMALL )
          warning("score is essentially zero in lattice\n");

	(*it).gmtk_score.setFromLogP(score);
      } while ( it.next() );
    }
  }

#if 0
  // print the lattice for debugging reasons
  printf("acscale=%f, amscale=%f, lmscale=%f\n", _acscale, _amscale, _lmscale);
  printf("starting %d and end %d\n", _start, _end);
  for ( unsigned i = 0; i < _numberOfNodes; i++ ) {
    printf("node %d at frame (%u,%u):\n", i, _latticeNodes[i].startFrame, _latticeNodes[i].endFrame);
    if ( _latticeNodes[i].edges.totalNumberEntries() > 0 ) {
      shash_map_iter<unsigned, LatticeEdge>::iterator it;
      _latticeNodes[i].edges.begin(it);
      do {
	LatticeEdge &edge = *it;
	printf("\tto %u, w=%d, ac=%f, lm=%f, p=%f, t=%f gmtk=%f\n", it.key(), edge.emissionId, edge.ac_score.val(), edge.lm_score.val(), edge.posterior.val(), _latticeNodes[it.key()].time, edge.gmtk_score.val());
      } while ( it.next() );
    }
  }
#endif
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::normalizePosterior
 *     in HTK lattices, poteriors are proporgated:
 *     For each node, the sum poteriors of outgoing links equal to the poterior
 *     of incoming edge.  This is not best suitable for GMTK because
 *     we need the sum posteriors of all outgoing links to be unity.
 *-----------------------------------------------------------------------
 */
void LatticeADT::normalizePosterior() {
  // normalizing posteriors has the benefits of better prunning.
  // We don't need to do anything on the gmtk_score because
  // those contains the real log probabilities which can be used
  // with proper scale factors.
  for ( unsigned i = 0; i < _numberOfNodes; ++i ) {
    if ( _latticeNodes[i].edges.totalNumberEntries() > 0 ) {
      shash_map_iter<unsigned, LatticeEdge>::iterator it;
      _latticeNodes[i].edges.begin(it);
      logpr sum;
      sum.set_to_zero();
      do {
	sum += (*it).posterior;
      } while ( it.next() );

      _latticeNodes[i].edges.begin(it);
      do {
	(*it).posterior = (*it).posterior / sum;
      } while ( it.next() );
    }
  }
}


const double LatticeADT::_frameRate = 100.0;
