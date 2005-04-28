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

#include "GMTK_LatticeADT.h"
#include "GMTK_NamedObject.h"
#include "GMTK_Vocab.h"
#include "GMTK_GMParms.h"
#include "error.h"


LatticeADT::LatticeADT() : _latticeNodes(NULL), _frameRate(100.0) {
}


LatticeADT::~LatticeADT() {
	delete [] _latticeNodes;
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

	_latticeNodes = new LatticeNode [_numberOfNodes];
	char *line = new char [1024];
	const char seps[] = " \t\n";
	unsigned id;

	// reading nodes
	unsigned frame;
	double time;
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
				time = atof(ptr+2);
				frame = (unsigned) (time * _frameRate);
				_latticeNodes[id].startFrame = _latticeNodes[id].endFrame = frame;
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
				break;
			default:
				error("unknonw token %s in LatticeCTP::readFromHTKLattice\n", ptr);
				break;
			}
		}

		_latticeNodes[id].edges.insert(endNodeId, edge);
	}

	delete [] line;
	
	// print the lattice
	for ( unsigned i = 0; i < _numberOfNodes; i++ ) {
		printf("node %d at frame (%d,%d):\n", i, _latticeNodes[i].startFrame, _latticeNodes[i].endFrame);
		if ( _latticeNodes[i].edges.totalNumberEntries() > 0 ) {
			shash_map2<unsigned, LatticeEdge>::iterator it = _latticeNodes[i].edges.begin();
			do {
				LatticeEdge &edge = *it;
				printf("\tto %d, w=%d, ac=%f, lm=%f\n", it.key(), edge.emissionId, edge.ac_score.val(), edge.lm_score.val());
			} while ( it.next() );
		}
	}
}


/*-
 *-----------------------------------------------------------------------
 * LatticeADT::read(ifs)
 *     read the lattice object from master file
 *
 * Results:
 *     no results:
 *-----------------------------------------------------------------------
 */
void LatticeADT::read(iDataStreamFile &is) {
	// read in name
	NamedObject::read(is);

	// read in node cardinality
	is.read(_nodeCardinality, "Can't read node cardinality");

	// read in lattice filename
	char *latticeFile;
	is.read(latticeFile, "Can't read lattice filename");

	// read in vocab filename
	string vocabName;
	is.read(vocabName, "Can't read vocab name");

	if ( GM_Parms.vocabsMap.find(vocabName) == GM_Parms.vocabsMap.end() )
		error("Error: reading file '%s' line '%d', LatticeCPT '%s' specifies Vocab name '%s' that does not exist", is.fileName(), is.lineNo(), _name.c_str(), vocabName.c_str());

	Vocab *vocab = GM_Parms.vocabs[GM_Parms.vocabsMap[vocabName]];
	iDataStreamFile lfifs(latticeFile);
	readFromHTKLattice(lfifs, *vocab);
	_wordCardinality = vocab->size();
}


#ifdef MAIN


#include "rand.h"
#include "GMTK_ObservationMatrix.h"
#include "GMTK_GMParms.h"


RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;
 

int main() {
	LatticeADT cpt;
	Vocab vocab(13);
	vocab.read("latticeVocab.txt");
	
	iDataStreamFile ifs("htk_sample.lat");
	cpt.readFromHTKLattice(ifs, vocab);

	return 0;
}


#endif

