//=========================================================================================

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "GMTK_Vocab.h"
#include "GMTK_FNGramCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_FileParser.h"
#include "GMTK_HidDiscRV.h"
#include "fileParser.h"
#include "logp.h"
#include "rand.h"
#include "error.h"


RAND rnd;
GMParms GM_Parms;
ObservationMatrix globalObservationMatrix;

//=========================================================================================

void testFLM1 ( char	*masterFileName,
		char	*strFileName ) 
{

//	fprintf ( stderr, " **** in testFLM1 **** \n" );
	// read in master file
	{
//		fprintf ( stderr, " **** read in input master file <%s> ... \n", masterFileName );
		iDataStreamFile pf(masterFileName, false, true);
		GM_Parms.read(pf);
	}
	GM_Parms.finalizeParameters();
	if ( GM_Parms.fngramImpsMap.find(string("oc_3tags")) == GM_Parms.fngramImpsMap.end() ) {
		fprintf ( stderr, " **** cannot find oc_3tags in FNGramImp\n");
		return;
	} else {
//		fprintf ( stderr, " --> found oc_3tags in FNGramImp ... \n" );
	}

	// read in structure file
//	fprintf ( stderr, " **** read in model structure file ... \n" );
	FileParser fp(strFileName);

	// parse the file
//	fprintf ( stderr, " **** parse graphical model ... \n" );
	fp.parseGraphicalModel();
	// create the rv variable objects
//	fprintf ( stderr, " **** create the RV objects ... \n" );
	fp.createRandomVariableGraph();
	// Make sure that there are no directed loops in the graph.
//	fprintf ( stderr, " **** ensure valid template ... \n" );
	fp.ensureValidTemplate();
//	fprintf ( stderr, " **** associate with data params ... \n" );
	fp.associateWithDataParams(FileParser::noAllocate);

	// get ngram pointer
//	fprintf ( stderr, " **** looking for fngrams ... \n" );
//	fprintf ( stderr, "number of fngram imps %d\n", GM_Parms.fngramImps.size());
//	fprintf ( stderr, "number of fngrams %d\n", GM_Parms.fngramCpts.size());
	if ( GM_Parms.fngramCptsMap.find(string("oc_3tags")) == GM_Parms.fngramCptsMap.end() ) {
		fprintf ( stderr, " **** cannot find oc_3tags in FNGramCPT\n");
		return;
	} else {
//		fprintf ( stderr, " --> found oc_3tags in FNGramCPT ... \n" );
	}

	FNGramCPT *fngram = GM_Parms.fngramCpts[GM_Parms.fngramCptsMap[string("oc_3tags")]];

//	fprintf ( stderr, "**** fngram pointer: %x\n", fngram);
	if ( fngram == NULL ) {
		fprintf ( stderr, " ERROR null pointer ... ??? \n" );
		return;
	}

	// read in vocabulary

//	fprintf ( stderr, " reading in tagVocab.list ... \n" );
	Vocab voc_tag(50);
	voc_tag.read("./data/tagVocab.list");

//	fprintf ( stderr, " reading in binaryB2.list ... \n" );
	Vocab voc_oc(4);
	voc_oc.read("./data/binaryB2.list");

//	fprintf ( stderr, " setting up RV information ... \n" );

	fprintf ( stderr, " PREPARING TO TEST ... \n" );

	RVInfo		ocInfo, tagInfo;
	HidDiscRV	oc2(ocInfo,2,4);
	HidDiscRV	tag0(tagInfo,2,50), tag1(tagInfo,1,50), tag2(tagInfo,0,50);

	vector<RV*>	parents(3);
	parents[0] = &tag2;
	parents[1] = &tag1;
	parents[2] = &tag0;

	logpr prob;

#if 1

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-comma, T-bos, T-bos ]
	tag2.val = voc_tag.index("T-comma");
	tag1.val = voc_tag.index("<s>");
	tag0.val = voc_tag.index("T-bos");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );

	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

//	exit (-1);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-wrb, T-wrb, T-nns ]
	tag2.val = voc_tag.index("T-wrb");
	tag1.val = voc_tag.index("T-wrb");
	tag0.val = voc_tag.index("T-nns");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-<s>, T-<s>, T-<s> ]
	tag2.val = voc_tag.index("T-<s>");
	tag1.val = voc_tag.index("T-<s>");
	tag0.val = voc_tag.index("T-<s>");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-period, T-ls, T-ls ]
	tag2.val = voc_tag.index("T-period");
	tag1.val = voc_tag.index("T-ls");
	tag0.val = voc_tag.index("T-ls");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-ls, T-period, T-ls ]
	tag2.val = voc_tag.index("T-ls");
	tag1.val = voc_tag.index("T-period");
	tag0.val = voc_tag.index("T-ls");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-ls, T-ls, T-period ]
	tag2.val = voc_tag.index("T-ls");
	tag1.val = voc_tag.index("T-ls");
	tag0.val = voc_tag.index("T-period");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-cc, T-bos, T-ls ]
	tag2.val = voc_tag.index("T-cc");
	tag1.val = voc_tag.index("T-bos");
	tag0.val = voc_tag.index("T-ls");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-cc, T-ls, T-bos ]
	tag2.val = voc_tag.index("T-cc");
	tag1.val = voc_tag.index("T-ls");
	tag0.val = voc_tag.index("T-bos");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-ls, T-cc, T-bos ]
	tag2.val = voc_tag.index("T-ls");
	tag1.val = voc_tag.index("T-cc");
	tag0.val = voc_tag.index("T-bos");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);

	// --------------------------------------------------------------------------------
	// Pr [ B-0 | T-cc, T-cc, T-bos ]
	tag2.val = voc_tag.index("T-cc");
	tag1.val = voc_tag.index("T-cc");
	tag0.val = voc_tag.index("T-bos");
	oc2.val  = voc_oc.index("B-0");

	fprintf ( stderr, " computing Pr [ %d | %d, %d, %d ] \n",
		oc2.val, tag2.val, tag1.val, tag0.val );
	
	prob = fngram->probGivenParents ( parents, &oc2 );

	fprintf ( stderr, " --> probability = %f = %f \n\n", prob.val(), prob.val()/M_LN10);
#endif



#if 0
	// --------------------------------------------------------------------------------
	int	ii, jj, kk;

	oc2.val = voc_oc.index("B-0");

	fprintf ( stderr, " \n\n\n" );
	fprintf ( stderr, " STARTING BIG TRIPLE LOOP ... \n\n\n" );

	for ( ii=2; ii<49; ii++ ) {
		tag2.val = ii;
		for ( jj=2; jj<49; jj++ ) {
			tag1.val = jj;
			for ( kk=2; kk<49; kk++ ) {
				tag0.val = kk;
				prob = fngram->probGivenParents ( parents, &oc2 );
				fprintf ( stderr, " Pr [ B-0 | %2d, %2d, %2d ] = %7.3f \n", ii, jj, kk, prob.val() );
			}
		}
	}
#endif

}

//=========================================================================================

int main ( int	 argc,
	   char	*argv[] ) 
{
	char	masterFileName[64];
	char	strFileName[64];

	ieeeFPsetup();
	set_new_handler(memory_error);

	if ( argc == 3 ) {
		strcpy ( masterFileName, argv[1] );
		strcpy ( strFileName,    argv[2] );
	} else {
		fprintf ( stderr, " parameters required : masterFileName, strFileName \n" );
	}

	testFLM1 ( masterFileName, strFileName );

}

//=========================================================================================
