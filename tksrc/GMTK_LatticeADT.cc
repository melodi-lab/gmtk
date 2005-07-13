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
LatticeADT::LatticeADT() : _latticeNodes(NULL), _numberOfNodes(0),
			   _numberOfLinks(0), _start(0), _end(0),
			   _lmscale(0), _wdpenalty(0), _acscale(0),
			   _amscale(0), _base(0), _frameRelax(0),
			   _latticeFile(NULL), _numLattices(0), _curNum(0),
			   _nodeCardinality(0), _wordCardinality(0) {
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
 *     read an HTK lattice
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
void LatticeADT::readFromHTKLattice(iDataStreamFile &ifs, const Vocab &vocab) {
	char *s_tmp, *ptr;
	do {
		ifs.read(s_tmp);
		if ( s_tmp == NULL )
			break;

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
		} else if ( strstr(s_tmp, "amscale") != NULL ) {
			// parse am scale
			ptr = strchr(s_tmp, '=');
			_amscale = atof(++ptr);
		} else if ( strstr(s_tmp, "start=") != NULL ) {
			// parse start node id
			ptr = strchr(s_tmp, '=');
			_start = (unsigned)atoi(++ptr);
		} else if ( strstr(s_tmp, "end=") != NULL ) {
			// parse start node id
			ptr = strchr(s_tmp, '=');
			_end = (unsigned)atoi(++ptr);
		} else if ( strstr(s_tmp, "N=") != NULL ) {
			// parse number of nodes
			ptr = strchr(s_tmp, '=');
			_numberOfNodes = (unsigned)atoi(++ptr);
			if ( _numberOfNodes > _nodeCardinality )
				error("Runtime error: number of nodes in lattice %d is bigger than node cardinlaity %d", _numberOfNodes, _nodeCardinality);
		} else if ( strstr(s_tmp, "L=") != NULL ) {
			// parse number of nodes
			ptr = strchr(s_tmp, '=');
			_numberOfLinks = (unsigned)atoi(++ptr);
			break;		// read to parse node and link information
		}
		
		delete [] s_tmp;
		s_tmp = NULL;
	} while ( ! ifs.isEOF() );

	if ( _latticeNodes )
		delete [] _latticeNodes;
	_latticeNodes = new LatticeNode [_numberOfNodes];
	char *line = new char [1024];
	const char seps[] = " \t\n";
	unsigned id;

	// reading nodes
	for ( unsigned i = 0; i < _numberOfNodes; i++ ) {
		ifs.readLine(line, 1024);
		ptr = strtok(line, seps);
		if ( ptr[0] != 'I' || ptr[1] != '=' )
			error("expecting I= in line %s at LatticeCPT::readFromHTKLattice", ptr);
		if ( (id = atoi(ptr+2)) > _numberOfNodes )
			error("node id %d is bigger than number of nodes in LatticeCPT::readFromHTKLattice", id, _numberOfNodes);
		if ( (ptr = strtok(NULL, seps)) != NULL ) {
			// there is time information
			if ( ptr[0] == 't' ) {
				// set up the time
				_latticeNodes[id].time = atof(ptr+2);
				// by default, there is no frame constrain
				_latticeNodes[id].startFrame = 0;
				_latticeNodes[id].endFrame = ~0;
			}
		}
	}

	// reading links
	double score;
	unsigned endNodeId = 0;
	for ( unsigned i = 0; i < _numberOfLinks; i++ ) {
		ifs.readLine(line, 1024);
		ptr = strtok(line, seps);
		if ( ptr[0] != 'J' || ptr[1] != '=' )
			error("expecting J= in line %s at LatticeCPT::readFromHTKLattice", ptr);
		if ( (id = atoi(ptr+2)) > _numberOfLinks )
			error("link id %d is bigger than number of nodes in LatticeCPT::readFromHTKLattice", id, _numberOfLinks);

		// starting node
		if ( (ptr = strtok(NULL, seps)) == NULL )
			error("expect starting node for link %d", id);
		if ( ptr[0] != 'S' || ptr[1] != '=' )
			error("expect starting node for link %d", id);
		if ( (id = atoi(ptr+2)) > _numberOfNodes )
			error("node id %d is bigger than number of nodes in LatticeCPT::readFromHTKLattice", id, _numberOfNodes);

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
					warning("Warning: word '%s' in lattice '%s' cannot be found in vocab", ptr+2, ifs.fileName());
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
				// probability
				score = atof(ptr+2);
				edge.prob_score.setFromP(score);
				break;
			default:
				error("unknonw token %s in LatticeCTP::readFromHTKLattice\n", ptr);
				break;
			}
		}

		_latticeNodes[id].edges.insert(endNodeId, edge);
	}

	delete [] line;
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
	is.read(_latticeFileName, "Cano't read lattice file or number of fetures");

	// if lattice file name is a string, then it is a list of lattices
	if ( ! strIsInt(_latticeFileName.c_str(), (int*)&_nodeCardinality) ) {
		if ( iterable() )
			error("ERROR: in lattices named '%s' in file '%s' line %d, can't have lattices defined recursively in files", name().c_str(),is.fileName(),is.lineNo());

		initializeIterableLattice(_latticeFileName);
	} else {
		// read in lattice filename
		char *latticeFile;
		is.read(latticeFile, "Can't read lattice filename");

		// read in vocab filename
		string vocabName;
		is.read(vocabName, "Can't read vocab name");

		if ( GM_Parms.vocabsMap.find(vocabName) == GM_Parms.vocabsMap.end() )
			error("Error: reading file '%s' line '%d', LatticeCPT '%s' specifies Vocab name '%s' that does not exist", is.fileName(), is.lineNo(), _name.c_str(), vocabName.c_str());

		// set vocabulary and read in HTK lattice
		Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
		iDataStreamFile lfifs(latticeFile);
		readFromHTKLattice(lfifs, *vocab);
		_wordCardinality = vocab->size();

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
void LatticeADT::seek(unsigned nmbr) {
	if (!iterable())
		error("ERROR: trying to seek in non-iterable lattice, '%s'\n",
	_curName.c_str() );

	// rewind the file
	_latticeFile->rewind();
	_latticeFile->read(_numLattices, "num lattices");
	if ( _numLattices <= nmbr )
		error("number of lattices (%d) in '%s' is less than segment number (%d)", _numLattices, _latticeFileName.c_str(), nmbr);

	// one way to do this is like decision trees
	// here since every lattice cpt has same number of tokens
	// an easy approached is used in this implementation
	string tmp;
	for ( unsigned i = 0; i < nmbr; i++ ) {
		for ( unsigned j = 0; j < 6; j++ )
			_latticeFile->read(tmp);
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
void LatticeADT::nextIterableLattice() {
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

	// read in lattice filename
	char *latticeFile;
	_latticeFile->read(latticeFile, "Can't read lattice filename");

	// read in vocab filename
	string vocabName;
	_latticeFile->read(vocabName, "Can't read vocab name");

	if ( GM_Parms.vocabsMap.find(vocabName) == GM_Parms.vocabsMap.end() )
		error("Error: reading file '%s' line '%d', LatticeCPT '%s' specifies Vocab name '%s' that does not exist", _latticeFile->fileName(), _latticeFile->lineNo(), _name.c_str(), vocabName.c_str());

	// set vocabulary and read in HTK lattice
	Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
	iDataStreamFile lfifs(latticeFile);
	readFromHTKLattice(lfifs, *vocab);
	_wordCardinality = vocab->size();

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

	const LatticeNode *const end_p = _latticeNodes + _numberOfNodes;
	unsigned frame;
	for ( LatticeNode *p = _latticeNodes; p != end_p; ++p ) {
		frame = (unsigned)roundf(warpRate * p->time);
		p->startFrame = (frame < _frameRelax) ? 0 : frame - _frameRelax;
		frame += _frameRelax;
		p->endFrame = (frame > lastFrameId) ? lastFrameId : frame;
	}

#if 1
	// print the lattice for debugging reasons
	for ( unsigned i = 0; i < _numberOfNodes; i++ ) {
		printf("node %d at frame (%u,%u):\n", i, _latticeNodes[i].startFrame, _latticeNodes[i].endFrame);
		if ( _latticeNodes[i].edges.totalNumberEntries() > 0 ) {
			shash_map_iter<unsigned, LatticeEdge>::iterator it;
			_latticeNodes[i].edges.begin(it);
			do {
				LatticeEdge &edge = *it;
				printf("\tto %u, w=%d, ac=%f, lm=%f, p=%f, t=%f\n", it.key(), edge.emissionId, edge.ac_score.val(), edge.lm_score.val(), edge.prob_score.val(), _latticeNodes[it.key()].time);
			} while ( it.next() );
		}
	}
#endif
}

