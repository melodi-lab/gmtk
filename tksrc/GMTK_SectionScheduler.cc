/*
 * GMTK_SectionScheduler.h
 *   Root class for the inference algorithms at the time series level.
 *
 * Written by Richard Rogers <rprogers@uw.edu>
 *
 * Copyright (C) 2015 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include <algorithm>
#include <float.h>

#include "error.h"

#include "GMTK_BoundaryTriangulate.h"
#include "GMTK_CountIterator.h"

#include "GMTK_SectionScheduler.h"
#include "GMTK_SectionInferenceAlgorithm.h"
#include "GMTK_ObsDiscRV.h"

// default names of the three sections for printing/debugging messages.
const char* SectionScheduler::P1_n = "P'";
const char* SectionScheduler::Co_n = "C'";
const char* SectionScheduler::E1_n = "E'";


// clear all memory on each new segment.
bool SectionScheduler::perSegmentClearCliqueValueCache = true;
// turns on or off VE separators, and determins the type to use PC, or PCG. 
unsigned SectionScheduler::useVESeparators = (VESEP_PC | VESEP_PCG);
// turn off VE separators by default.
unsigned SectionScheduler::veSeparatorWhere = 0;
// uncomment to turn on VE separators.
// unsigned SectionScheduler::veSeparatorWhere = (VESEP_WHERE_P | VESEP_WHERE_C | VESEP_WHERE_E);
bool SectionScheduler::jtWeightUpperBound = false;
bool SectionScheduler::jtWeightMoreConservative = false;
float SectionScheduler::jtWeightPenalizeUnassignedIterated = 0.0;
float SectionScheduler::jtWeightSparseNodeSepScale = 1.0;
float SectionScheduler::jtWeightDenseNodeSepScale = 1.0;
const char* SectionScheduler::junctionTreeMSTpriorityStr = "DSU";
const char* SectionScheduler::interfaceCliquePriorityStr = "W";

bool SectionScheduler::viterbiScore = false;
bool SectionScheduler::onlineViterbi = false;
bool SectionScheduler::mmapViterbi = true;
bool SectionScheduler::sectionDoDist = false;
bool SectionScheduler::binaryViterbiSwap = false;

char * SectionScheduler::pVitTrigger = NULL;
char * SectionScheduler::cVitTrigger = NULL;
char * SectionScheduler::eVitTrigger = NULL;
bool   SectionScheduler::vitRunLength = false;

FILE * SectionScheduler::binaryViterbiFile = NULL;
char * SectionScheduler::binaryViterbiFilename = NULL;
off_t  SectionScheduler::binaryViterbiOffset;
off_t  SectionScheduler::nextViterbiOffset;
bool SectionScheduler::normalizePrintedCliques = true;







struct PairUnsigned1stElementCompare {  
  // sort ascending, comparing only 1st element.
  bool operator() (const pair<unsigned,unsigned>& a, 
		   const pair<unsigned,unsigned>& b) {
    // printf("comparing %d with %d\n",a.first,b.first);
    return (a.first) < (b.first);
  }
};




string
cliqueName(unsigned clique_num) {
  char buf[16];
  sprintf(buf, "C%03u", clique_num);
  return string(buf);
}

void
writeMsgOrder(oDataStreamFile &os, vector< pair<unsigned,unsigned> > const &msgs) {
  os.write(msgs.size()); os.nl();
  for (unsigned i=0; i < msgs.size(); ++i) {
    os.write(i); 
    os.write(cliqueName(msgs[i].first)); 
    os.write(cliqueName(msgs[i].second));     
    os.nl();
  }
}


void 
SectionScheduler::printAllIAInfo(char const* fileName, bool writeComments) {
  assert(fileName);

  oDataStreamFile os(fileName);
  os.setWriteCommentsStatus(writeComments);
  fp.writeGMId(os);
  gm_template.writeSections(os);

  os.nl();
  os.writeComment("---\n");
  os.writeComment("---- P inference architecture\n");
  P1.writeInferenceArchitecture(os, 'P');
  os.writeComment("---- within-section message order\n");
  writeMsgOrder(os, P1_message_order);
  os.writeComment("---- P has no left interface\n");
  os.write(0); os.nl();
  os.writeComment("---- P right interface\n");
  if (P1.cliques.size() > 0) {
    os.write(P_ri_to_C.size()); os.nl();
    for (unsigned i=0; i < P_ri_to_C.size(); ++i) {
      os.write(i); os.write(cliqueName(P_ri_to_C[i])); os.nl();
    }
  } else {
    os.write(0); os.nl();
  }
  os.nl();

  os.writeComment("---\n");
  os.writeComment("---- C inference architecture\n");
  Co.writeInferenceArchitecture(os, 'C');
  os.writeComment("---- within-section message order\n");
  writeMsgOrder(os, Co_message_order);
  os.writeComment("---- C left interface\n");
  os.write(C_li_to_C.size()); os.nl();
  for (unsigned i=0; i < C_li_to_C.size(); ++i) {
    os.write(i); os.write(cliqueName(C_li_to_C[i])); os.nl();
  }
  os.writeComment("---- C right interface\n");
  os.write(C_ri_to_C.size()); os.nl();
  for (unsigned i=0; i < C_ri_to_C.size(); ++i) {
    os.write(i); os.write(cliqueName(C_ri_to_C[i])); os.nl();
  }
  os.nl();

  os.writeComment("---\n");
  os.writeComment("---- E inference architecture\n");
  E1.writeInferenceArchitecture(os, 'E');
  os.writeComment("---- within-section message order\n");
  writeMsgOrder(os, E1_message_order);
  os.writeComment("---- E left interface\n");
  if (E1.cliques.size() > 0) {
    os.write(E_li_to_C.size()); os.nl();
    for (unsigned i=0; i < E_li_to_C.size(); ++i) {
      os.write(i); os.write(cliqueName(E_li_to_C[i])); os.nl();
    }
  } else {
    os.write(0); os.nl();
  }
  os.writeComment("---- E has no right interface\n");
  os.write(0); os.nl();
}



// Formerly JunctionTree::printAllJTInfo()
void 
SectionScheduler::printInferencePlanSummary(char const *fileName) {
  FILE* f;
  if ((fileName == NULL)
      ||
      ((f = ::fopen(fileName,"w")) == NULL))
    return;


  // print partition (clique,separator) information

  fprintf(f,"===============================\n");
  fprintf(f,"   P1 partition information: JT_weight = %f\n",
	  0.0 /*junctionTreeWeight(P1,P_ri_to_C,NULL,&Co.nodes)*/);
  printAllJTInfo(f,P1,P_ri_to_C,NULL,&Co.nodes);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");
  fprintf(f,"   Co partition information: JT_weight = %f\n",
	  0.0 /*junctionTreeWeight(Co,C_ri_to_E,&P1.nodes,&E1.nodes)*/);
  printAllJTInfo(f,Co,C_ri_to_C,&P1.nodes,&E1.nodes);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");
  fprintf(f,"   E1 partition information: JT_weight = %f\n",
	  0.0 /*junctionTreeWeight(E1,E_root_clique, &Co.nodes,NULL)*/);
  printAllJTInfo(f,E1,E_root_clique,&Co.nodes,NULL);
  fprintf(f,"\n\n");

  // print message order information
  fprintf(f,"===============================\n\n");    

  fprintf(f,"===============================\n");  
  fprintf(f,"   P1 message order\n");
  printMessageOrder(f,P1_message_order);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");  
  fprintf(f,"   Co message order\n");
  printMessageOrder(f,Co_message_order);
  fprintf(f,"\n\n");

  fprintf(f,"===============================\n");  
  fprintf(f,"   E1 message order\n");
  printMessageOrder(f,E1_message_order);
  fprintf(f,"\n\n");


  fclose(f);
}


void
SectionScheduler::printAllJTInfo(FILE* f,
			     JT_Partition& part,
			     vector<unsigned> const &roots,
			     set <RV*>* lp_nodes,
			     set <RV*>* rp_nodes)

{
  // print cliques information
  fprintf(f,"=== Clique Information ===\n");
  fprintf(f,"Number of cliques = %ld\n",(unsigned long)part.cliques.size());
  if (part.cliques.size() > 0)
    for (unsigned i=0; i < roots.size(); ++i)
      printAllJTInfoCliques(f,part,roots[i],0,lp_nodes,rp_nodes);

  // print separator information
  fprintf(f,"\n=== Separator Information ===\n");
  fprintf(f,"Number of separators = %ld\n",(unsigned long)part.separators.size());
  for (unsigned sepNo=0;sepNo<part.separators.size();sepNo++) {
    fprintf(f,"== Separator number: %d\n",sepNo);
    part.separators[sepNo].printAllJTInfo(f);
  }

}


void
SectionScheduler::printAllJTInfoCliques(FILE* f,
				    JT_Partition& part,
				    unsigned root,
				    const unsigned treeLevel,
				    set <RV*>* lp_nodes,
				    set <RV*>* rp_nodes)
{
  const bool is_P_partition = (lp_nodes == NULL);
  const bool is_E_partition = (rp_nodes == NULL);
  // can't be both.
  assert ( !(is_P_partition && is_E_partition) );

  // print clique's information
  for (unsigned i=0;i<treeLevel;i++) fprintf(f,"  ");
  fprintf(f,"== Clique number: %d",root);
  if (treeLevel == 0)
    fprintf(f,", root/right-interface clique");
  if (treeLevel == 0) {
    if (!is_E_partition)
      fprintf(f,", root/right-interface clique");
    else
      fprintf(f,", root clique");
  }
  if (part.cliques[root].ceReceiveSeparators.size() >
      part.cliques[root].numVESeparators() + part.cliques[root].children.size() ) {
    if (part.cliques[root].ceReceiveSeparators.size() > 1 + part.cliques[root].numVESeparators()) 
      fprintf(f,", left-interface clique");
    else
      fprintf(f,", leaf/left-interface clique");
  } else {
    assert ( part.cliques[root].ceReceiveSeparators.size() == 
	     part.cliques[root].children.size() + part.cliques[root].numVESeparators());
    if (part.cliques[root].ceReceiveSeparators.size() == part.cliques[root].numVESeparators()) 
      fprintf(f,", leaf");
  }
  fprintf(f,"\n");
  part.cliques[root].printAllJTInfo(f,treeLevel,part.unassignedInPartition,
				    jtWeightUpperBound,
				    jtWeightMoreConservative,
				    true,
				    lp_nodes,rp_nodes);
  for (unsigned childNo=0;
       childNo<part.cliques[root].children.size();childNo++) {
    unsigned child = part.cliques[root].children[childNo];
    printAllJTInfoCliques(f,part,child,treeLevel+1,lp_nodes,rp_nodes);
  }
}



void
SectionScheduler::printMessageOrder(FILE *f,
				vector< pair<unsigned,unsigned> >& message_order)
{
  fprintf(f,"Number of messages: %ld\n",(unsigned long)message_order.size());
  for (unsigned m=0;m<message_order.size();m++) {
    const unsigned from = message_order[m].first;
    const unsigned to = message_order[m].second;
    fprintf(f,"  %d: %d --> %d\n",m,from,to);
  }
}

  // Formerly GMTemplate::reportScoreStats()
void 
SectionScheduler::reportScoreStats() {
}


  // NOTE:  This assumes all inference algorithms will have clique-like things they
  //        need to print posteriors of

  // Set the range of selected clique #'s in P', C', E' for printing.
  // TODO: preconditions
void 
SectionScheduler::setCliquePrintRanges(char *p_range, char *c_range, char *e_range) {
  if (p_range != NULL) {
    BP_Range* tmp = new BP_Range(p_range,0,gm_template.P.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      p_clique_print_range = tmp;
  }
  if (c_range != NULL) {
    BP_Range* tmp = new BP_Range(c_range,0,gm_template.C.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      c_clique_print_range = tmp;
  }
  if (e_range != NULL) {
    BP_Range* tmp = new BP_Range(e_range,0,gm_template.E.cliques.size());
    if (!tmp || tmp->length() <= 0)
      delete tmp;
    else
      e_clique_print_range = tmp;
  }
}


  // Print to f the order of the variables in each clique selected by setCliquePrintRanges().
void 
SectionScheduler::printCliqueOrders(FILE *f) {
}

  // Returns the size (in # of floats) of the cliques selected by setCliquePrintRanges().
void 
SectionScheduler::getCliquePosteriorSize(unsigned &p_size, unsigned &c_size, unsigned &e_size) {
}


void 
SectionScheduler::printAllCliques(FILE *f,const bool normalize, const bool unlog,
				  const bool justPrintEntropy,
				  ObservationFile *obs_file)
{
  // FIXME - Why 2 iterators? in the old (1.X) GMTK, setCurrentInferenceShiftTo()
  //         adjusted the JunctionTree::inference_it member. The SectionScheduler
  //         no longer has that member, so setCurrentInferenceShiftTo() takes the
  //         SectionIterator to adjust as an argument. But in this code, the current
  //         section is incremented via stss_it++ operations, thus 
  //         setCurrentInferenceShiftTo(stss_it, stss_it.cur_st()) would find the
  //         iterator already at the desired section and therefore skip updating the
  //         RV frames. As a temporary hack, I added the second inference_it so it
  //         can train a section behind stss_it and thus trigger the update. It would
  //         be better to just iterate s/currentPartition/current_section/ over the
  //         section_table_array indices.
  SectionIterator stss_it(*this), inference_it(*this);

  if (p_clique_print_range != NULL) {
    setCurrentInferenceShiftTo(inference_it, stss_it.cur_st());
    section_table_array[stss_it.cur_st()]->printAllCliques(f, p_clique_print_range, stss_it, 
							   section_structure_array[stss_it.cur_ss()],
							   normalize, unlog, justPrintEntropy, obs_file);
  }
  stss_it++;

  if (c_clique_print_range != NULL) {
    for (; !stss_it.at_last_entry(); stss_it++) {
      int currentPartition = stss_it.cur_st();
      setCurrentInferenceShiftTo(inference_it, currentPartition);
      section_table_array[stss_it.cur_st()]->printAllCliques(f, c_clique_print_range, stss_it,
							     section_structure_array[stss_it.cur_ss()],
							     normalize, unlog, justPrintEntropy, obs_file);
    }
  } else {
    // partNo = section_structure_array.size()-1;
    stss_it.set_to_last_entry();
  }

  if (e_clique_print_range != NULL) {
    setCurrentInferenceShiftTo(inference_it, stss_it.cur_st());
    section_table_array[stss_it.cur_st()]->printAllCliques(f, e_clique_print_range, stss_it,  
							   section_structure_array[stss_it.cur_ss()],
							   normalize, unlog, justPrintEntropy, obs_file);
  }
}



/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::prepareForUnrolling()
 *
 *  Sets up a number of data structures that must be ready before
 *  we can do general unrolling. This includes:
 *  
 *   1) computes the number of hidden variables in each clique/separator
 *      (i.e., the number of discrete hidden and non-continous variables)
 *
 * Preconditions:
 *   computeSeparatorIterationOrders() must have been called.
 *
 * Postconditions:
 *   everything is now ready for general unrolling.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques/separators within each section.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::prepareForUnrolling() {
  const unsigned totalNumVESeps =
    P1.numVEseps + Co.numVEseps + E1.numVEseps;
  if (useVESeparators && totalNumVESeps > 0) {
    // possibly set VE sep file

    SeparatorClique::veSeparatorFile = NULL;

    if (SeparatorClique::recomputeVESeparatorTables || 
	((SeparatorClique::veSeparatorFileName != NULL) &&
	 (::fsize(SeparatorClique::veSeparatorFileName) == 0))) {
      // open file for writing since we're re-generating the information.

      if (SeparatorClique::veSeparatorFileName != NULL) {
	SeparatorClique::veSeparatorFile =
	  ::fopen(SeparatorClique::veSeparatorFileName,"w");
	if (SeparatorClique::veSeparatorFile == NULL) {
	  error("ERROR: cannot open VE separator file (%s) for writing.",SeparatorClique::veSeparatorFileName);
	}
      }
      SeparatorClique::generatingVESeparatorTables = true;
      infoMsg(IM::Default,"Computing information for %d total VE separators\n",totalNumVESeps);
    } else if (SeparatorClique::veSeparatorFileName != NULL) {
      // assume that the current ve sep file is valid.
      SeparatorClique::veSeparatorFile =
	::fopen(SeparatorClique::veSeparatorFileName,"r");
      if (SeparatorClique::veSeparatorFile == NULL) {
	error("ERROR: cannot open VE separator file (%s) for reading.",SeparatorClique::veSeparatorFileName);
      }
      SeparatorClique::generatingVESeparatorTables = false;
      infoMsg(IM::Default,"Reading information for %d total VE separators\n",totalNumVESeps);
    }
  }

  if (P1.numVEseps > 0) {
    infoMsg(IM::Default,"Computing P's %d VE separators\n",P1.numVEseps);
  }
  prepareForUnrolling(P1);

  if (Co.numVEseps > 0) {
    infoMsg(IM::Default,"Computing C's %d VE separators\n",Co.numVEseps);
  }
  prepareForUnrolling(Co);

  if (E1.numVEseps > 0) {
    infoMsg(IM::Default,"Computing E's %d VE separators\n",E1.numVEseps);
  }
  prepareForUnrolling(E1);

  if (useVESeparators && totalNumVESeps > 0) {
    // close file
    if (SeparatorClique::veSeparatorFile != NULL)
      fclose(SeparatorClique::veSeparatorFile);
    if (SeparatorClique::generatingVESeparatorTables)
      infoMsg(IM::Default,"Done computing information for %d total VE separators\n",totalNumVESeps);
  }
}

void 
SectionScheduler::prepareForUnrolling(JT_Partition &section) {
  for (unsigned i=0;i<section.cliques.size();i++) {
    section.cliques[i].prepareForUnrolling();
  }
  for (unsigned i=0;i<section.separators.size();i++) {
    section.separators[i].prepareForUnrolling();
  }
}


void checkDisps(sArray< MaxClique::AssignedNodeDisposition > const &d) {
  for (int i=0; i < d.len(); ++i) {
    assert(0 <= d[i] && d[i] <= 5);
  }
}

void checkDisps(Section const &s) {
  for (unsigned i=0; i < s.cliques.size(); ++i) {
    checkDisps(s.cliques[i].dispositionSortedAssignedNodes);
  }
}

  // Initialize stuff at the model-level. See prepareForSegment() for segment-level initialization.
  // TODO: explain parameters
void 
SectionScheduler::setUpDataStructures(iDataStreamFile &tri_file,
				      char const *varSectionAssignmentPrior,
				      char const *varCliqueAssignmentPrior,
				      bool checkTriFileCards)
{

  // Utilize both the section information and elimination order
  // information already computed and contained in the file. This
  // enables the program to use external triangulation programs,
  // where this program ensures that the result is triangulated
  // and where it reports the quality of the triangulation.

  infoMsg(IM::Max,"Reading triangulation file '%s' ...\n", tri_file.fileName());
  if (!fp.readAndVerifyGMId(tri_file, checkTriFileCards)) {
    error("ERROR: triangulation file '%s' does not match graph given in structure file '%s'\n",
	  tri_file.fileName(),fp.fileNameParsing.c_str());
  }
  gm_template.readPartitions(tri_file);
  gm_template.readInferenceArchitectures(tri_file);

  infoMsg(IM::Max,"Triangulating graph...\n");

#if 0
  gm_template.triangulateSectionsByInferenceArchitecture();
#else
  // TODO: It looks like this is really just adding the edges to make the
  //       cliques specified in the tri file actual cliques in the "graph"
  //       by adding any missing edges. 
  gm_template.triangulatePartitionsByCliqueCompletion();
#endif

  if (1) { 
    // check that graph is indeed triangulated.
    // TODO: perhaps take this check out so that inference code does
    // not need to link to the triangulation code (either that, or put
    // the triangulation check in a different file, so that we only
    // link to tri check code).

    // TODO: Post-refactor, we're not assuming the graph must be triangulated?
    //       Move this into the subset of section inference algorithms that require it.
    BoundaryTriangulate triangulator(fp,
				     gm_template.maxNumChunksInBoundary(),
				     gm_template.chunkSkip(),1.0);
    triangulator.ensurePartitionsAreChordal(gm_template);
  }

  ////////////////////////////////////////////////////////////////////
  // CREATE JUNCTION TREE DATA STRUCTURES
  infoMsg(IM::Default,"Creating Junction Tree\n"); fflush(stdout);
  setUpJTDataStructures(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  prepareForUnrolling();
  infoMsg(IM::Default,"DONE creating Junction Tree\n"); fflush(stdout);
  ////////////////////////////////////////////////////////////////////
}



// gmtkMFA version
void 
SectionScheduler::setUpDataStructures(char const *varSectionAssignmentPrior,
				      char const *varCliqueAssignmentPrior)
{
  createSectionJunctionTrees(junctionTreeMSTpriorityStr);
#if 1
  computeSectionInterfacesMFA();
  createFactorCliques();
  createDirectedGraphOfCliques();

  assignRVsToCliques(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  assignFactorsToCliques();
  // TODO: assignScoringFactorsToCliques();
  setUpMessagePassingOrders();
  // create seps and VE seps.
  createSeparators();
  computeSeparatorIterationOrders();

  // -- -- used only to compute weight.
  getCumulativeUnassignedIteratedNodes(); 
  // -- --

  sortCliqueAssignedNodesAndComputeDispositions(varCliqueAssignmentPrior);
#endif
}

/*
 *  init_CC_CE_rvs(it)
 *    given an iterator representing an unrolled segment,
 *    we create the pair of random variables CC and CE that
 *    will be shifted around in time/position in order
 *    to do inference.
 *
 * Preconditions:
 *    1) section_structure_array must be set up and initialized before
 *    this routine is called.
 *    2) All random variables should currently be shifted to their 
 *       zero position (otherwise thigns will be totally out of sync).
 *
 */
void
SectionScheduler::init_CC_CE_rvs(SectionIterator &it)
{
  if (it.num_c_sections() > 2) { 
    // we only need to fill this set of we have more than 2 table
    // sections (meaning we will need to do some RV adjustment).
    set <RV*> tmp1,tmp2;
    tmp1 = section_structure_array[1].returnRVsAndTheirObservedParentsAsSet();
    tmp2 = section_structure_array[2].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CC_rvs);
    tmp1 = section_structure_array[3].returnRVsAndTheirObservedParentsAsSet();
    unionRVs(tmp1,tmp2,cur_CE_rvs);
  } else {
    cur_CC_rvs.clear();
    cur_CE_rvs.clear();
  } 
  cur_cc_shift = cur_ce_shift = 0;
  it.go_to_section_no(0);
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables CC so that the right position (the 2nd
 *  C of CC) is now at position 'pos'. If the CE is shifted,
 *  we first shift CE back to position 0 before shifting CC.
 */
void
SectionScheduler::shiftCCtoPosition(int pos, bool reset_observed)
{
  if (cur_ce_shift != 0) {
    assert ( cur_cc_shift == 0 );
    // we need to get CE back to zero.
    adjustFramesBy (cur_CE_rvs, -cur_ce_shift * gm_template.chunkSkip()*fp.numFramesInC(),
		    reset_observed);  // gmtkOnline can't observe old values - ticket #468
    cur_ce_shift = 0;
  }
  int delta = pos - cur_cc_shift;
  adjustFramesBy (cur_CC_rvs, delta * gm_template.chunkSkip()*fp.numFramesInC());
  cur_cc_shift = pos;
}

/*
 * shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
 *  Shift the random variables in 'rvs' in time (frames) by
 *  the difference between 'prevPos' and 'pos'. This is used  in 
 *  printSavedViterbiValues to adjust the rvs' frame numbers in
 *  the unpacking buffers.
 */

void
SectionScheduler::shiftOriginalVarstoPosition(vector<RV*> rvs, int pos, int &prevPos)
{
  set<RV*> uprvs(rvs.begin(),rvs.end());
  int delta = (pos - prevPos);
  adjustFramesBy(uprvs, delta);
  prevPos = pos;
}

/*
 * shiftCCtoPosition(int pos)
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables CE so that the right position (the E
 *  of CE) is now at position 'pos'. If the CC pair is shifted,
 *  we first shift CC back to position 0 before shifting CE.
 */
void
SectionScheduler::shiftCEtoPosition(int pos, bool reset_observed)
{
  if (cur_cc_shift != 0) {
    assert ( cur_ce_shift == 0 );
    // we need to get CC back to zero.
    adjustFramesBy (cur_CC_rvs, -cur_cc_shift * gm_template.chunkSkip()*fp.numFramesInC(),
                    reset_observed);  // gmtkOnline can't observe old values - ticket #468
    // FIXME - maybe reset_observed && !onlineViterbi
    cur_cc_shift = 0;
  }
  int delta = pos - cur_ce_shift;
  adjustFramesBy (cur_CE_rvs, delta * gm_template.chunkSkip()*fp.numFramesInC());
  cur_ce_shift = pos;
}



/*
 * setCurrentInferenceShiftTo(int pos)
 *
 *  This activates sections 'pos-1' and 'pos' so the
 *  right of two successive sections is active at position 'pos'.
 *  This routine could be called:
 *
 *      activateTwoAdjacentSectionsWithRightSectionAtPos(pos)
 *      alignRightSideOfSectionPairToPos(pos)
 * 
 *  Given an absolute section number 'pos', we shift the pair of
 *  random variables (either CC or CE) so the *right* section (the 2nd
 *  C of CC or the E of CE) is now at position 'pos'. This means that
 *  messages entirely within section 'pos' and entirely within
 *  section 'pos-1' and between sections 'pos-1' and 'pos' will be
 *  correct, but no other messages are guaranteed to be correct.
 *
 * preconditions:
 *    the SectionIterator must be initilzed for the
 *    current and appropriate segment length.
 *
 * side effects
 *   - modifies the random variables corresponding to a section pair
 *   - modifies the section iterator.
 *
 */
void
SectionScheduler::setCurrentInferenceShiftTo(SectionIterator &it, int pos)
{

  if (it.at_entry(pos)) {
    //    printf("== Already at shift %d\n",pos);
    return;
  }

  //   printf("========================================\n");
  //   printf("== Setting current inference shift to %d\n",pos);
  //   printf("========================================\n");

  it.go_to_section_no(pos);
  if (it.num_c_sections() <= 2) {
    // then do nothing, since nothing is never shifted away from zero.
  } else {
    // need to do some work. We need to get the right
    // of the appropriate pair (either the second C of CC in PCCE,
    // or the E of CE in PCCE) to be at position pos.
    if (it.at_p()) {
      // P has an intersection with the first C in PCCE, so we need
      // to get that first C into the right position. We do that
      // by getting both CC's into the right position.
      shiftCCtoPosition(0);
    } else if (it.at_e()) {
      // get the CE into the right position
      shiftCEtoPosition(pos - 3);
    } else {
      // Get CC into the right position, where 'pos' corresponds
      // to the second C in PCCE.
      if (pos <= 2)
	shiftCCtoPosition(0);
      else 
	shiftCCtoPosition(pos-2);
    }
  }
}









/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::unroll()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k times, k >= 0. (i.e., we have
 *   a graph that has (k+1) copies of C') and place them into this
 *   class's members section_structure_array and section_table_array.
 *
 *   This also set up data structures for inference. Some of what this routine
 *   does is copy data out of STL structures into faster customized structures
 *   (such as hash tables, arrays, etc.).
 *
 *   Sets things up so that inference can take place.
 *   Unrolling is in terms of new C' sections (so that if S=1, k corresponds to frames)
 *   but in general the network will be T(P') + k*T(C') + T(E) frames long, and in general
 *   T(C') >= S. This also depends on if the boundary algorithm was run or not. 
 *
 * Preconditions:
 *   The sections must be validly instantiated with cliques, and
 *   the routine assignRVsToCliques() must have been called.
 *
 * Postconditions:
 *   The sections are such that they now know who they are connected to on
 *   the right and left.
 *   Also, unrolled data structures are set up and ready for inference for this length.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   section.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
unsigned
SectionScheduler::unroll(const unsigned int numFrames,
			 const UnrollTableOptions tableOption,
			 unsigned *totalNumberSections)
{
  // note: the argument name is numFrames, which indicates
  // the number of frames in an observation file segment.


  ////////////////////////////////////////////////////////////////////////
  // various other parameters regarding the current segment. These
  // are all computed by SectionScheduler::unroll(...)
  // 
  // 1) The number of frames given, but in units of the amount by which
  // the basic template would need to be unrolled to match the given
  // number of frames.
  unsigned basicTemplateMaxUnrollAmount;
  // 2) For the number of frames given, the min amount by which we can
  // unroll and still re-use the rvs accross time.
  unsigned basicTemplateMinUnrollAmount;
  // 3) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given. This
  // value could be negative.
  int modifiedTemplateMaxUnrollAmount;
  // 4) how much the modified [P' C' E'] template would need to be
  // unrolled to match the number of frames given for the structure
  // case (i.e., this would either be [P' E'], [P' C' E'], or [P' C' C' E']
  // so this value is never more than 1 and could be negative.
  int modifiedTemplateMinUnrollAmount;
  // 5) How many usable frames are there in the current segment (for
  //    some boundaries/triangulations and their resulting sections,
  //    it might not be possible to have a segment that is a multiple
  //    of one (1) and so this number is guaranteed to be this.
  unsigned numUsableFrames;
  // 6) The offset by which we should start the observations (should be given
  // to observation matrix). Why this exists, see 5).
  unsigned frameStart;

#if 0
  // FIXME  delete this if unneeded

  // 7) The current total number of real frames.
  unsigned curNumFrames;
#endif
  ////////////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////////////
  // collect/distribute/probE support variables used in various routines.
  // current set of unrolled random variables
  vector <RV*> cur_unrolled_rvs;
  // current mapping from 'name+frame' to integer index into unrolled_rvs.
  map < RVInfo::rvParent, unsigned > cur_ppf;
  ////////////////////////////////////////////////////////////////////////

  unsigned M = gm_template.maxNumChunksInBoundary();
  unsigned S = gm_template.chunkSkip();

  // first create the unrolled set of random variables corresponding to this JT.
  if (!gm_template.computeUnrollParameters(numFrames,
					   basicTemplateMaxUnrollAmount,
					   basicTemplateMinUnrollAmount,
					   modifiedTemplateMaxUnrollAmount,
					   modifiedTemplateMinUnrollAmount,
					   numUsableFrames,
					   frameStart))
    error("Segment of %d frames too short with current GMTK template of length [P=%d,C=%d,E=%d] %d frames, and M=%d,S=%d boundary parameters. Use segments of at least P+M*C+E = %d frames, different template, or decrease M,S if >1. You can identify which segment(s) in the input file(s) are too short with a command like \"obs-info -p -i1 file1 ... | awk '$2 < %d {print $1}'\"\n",
	  numFrames,
	  fp.numFramesInP(),fp.numFramesInC(),fp.numFramesInE(),
	  fp.numFrames(),
	  M,S, 
	  fp.numFramesInP() + M * fp.numFramesInC() + fp.numFramesInE(),
	  fp.numFramesInP() + M * fp.numFramesInC() + fp.numFramesInE());
  // unused const int numCoSectionTables = modifiedTemplateMaxUnrollAmount+1;
  const int numCoSectionStructures= modifiedTemplateMinUnrollAmount+1;
  const int numSectionStructures= modifiedTemplateMinUnrollAmount+3;

  // we should never have more than 2 C structures since 2 is the max
  // necessary.
  assert ( numCoSectionStructures <= 2 );

  infoMsg(IM::Info,"Number of Current Frames = %d, Number of Currently Usable Frames = %d\n",numFrames,numUsableFrames);
  infoMsg(IM::Tiny,"Number Of Frames = %d, Unrolling Basic Template %d times, Modified Template %d times\n",
	  numFrames,
	  basicTemplateMaxUnrollAmount,
	  modifiedTemplateMaxUnrollAmount);

  // unrolled random variables
  // vector <RV*> unrolled_rvs;
  // mapping from 'name+frame' to integer index into unrolled_rvs.
  // map < RVInfo::rvParent, unsigned > ppf;

  // Only unroll the minimum amount since unrolling is expensive and
  // takes up much memory.
  fp.unroll(basicTemplateMinUnrollAmount,cur_unrolled_rvs,cur_ppf);

  // set the observed variables for now, but these may/will be modified later.
  setObservedRVs(cur_unrolled_rvs);

  prepareForNextInferenceRound();
  // this clears the shared caches in the origin cliques
  clearCliqueSepValueCache(perSegmentClearCliqueValueCache);

  // clear out the old and pre-allocate for new size.
  section_structure_array.clear();
  section_table_array.clear();

  // 
  // For the structure array, we always only need to unroll by the
  // amount corresponding to the amount by which the basic random
  // variables have been unrolled.
  section_structure_array.resize(numSectionStructures);
  unsigned sectionNo = 0;
  new (&section_structure_array[sectionNo++]) PartitionStructures(P1
							       ,cur_unrolled_rvs
							       ,cur_ppf
							       ,0*S*fp.numFramesInC(),
							       false // does not have a li separator
							       );
  for (int p=0;p<numCoSectionStructures;p++) {
    new (&section_structure_array[sectionNo]) PartitionStructures(Co,cur_unrolled_rvs,cur_ppf,p*S*fp.numFramesInC());
    sectionNo++;
  }
  new (&section_structure_array[sectionNo++]) 
    PartitionStructures(E1,cur_unrolled_rvs,cur_ppf,
			modifiedTemplateMinUnrollAmount*S*fp.numFramesInC());

   
  // 
  // re-allocate. We add three to modifiedTemplateUnrollAmount
  // since: 
  //  1) when we unroll by n we get n+1 copies (+1) and
  //     modifiedTemplateUnrollAmount is how much to unroll,
  //     but below is how much to allocate.
  //  2) we need a section for P, even if it is empty. (+1) 
  //  3) we need a section for E, even if it is empty. (+1)
  if (tableOption == LongTable) {
#if 0
    // then we pre-allocate the table array to correspond to the
    // entire length of the segment. I.e., each table is unique.  With
    // this option, it is possible to store all of the tables in
    // memory at the same time (assuming there is sufficient RAM).
    section_table_array.resize(modifiedTemplateMaxUnrollAmount+3);
    unsigned sectionNo = 0;
    new (&section_table_array[sectionNo++]) PartitionTables(P1);
    for (int p=0;p<numCoSectionTables;p++) {
      new (&section_table_array[sectionNo]) PartitionTables(Co);
      sectionNo++;
    }
    new (&section_table_array[sectionNo++])
      PartitionTables(E1);
    assert (sectionNo == section_table_array.size());
#else

    // s/PartitionTables/SectionTables/ is now a class hierarchy to 
    // permit different SectionInferenceAlgorithm subclasses to do
    // inference in each section. So, we can't instantiate the actual
    // section_table_array entries until the particular SectionInferenceAlgorithm
    // for the section in question has been chosen.

    section_table_array.resize(modifiedTemplateMaxUnrollAmount+3);
    for (unsigned p=0; p < section_table_array.size(); ++p) {
      section_table_array[p] = NULL;
    }
#endif
  } else  if (tableOption == ShortTable) {
#if 0
    // Then we pre-allocate the table array to correspond only to the
    // structure array, assuming that the tables are going to be
    // reused in some way. In this option, it is not possible to store
    // all the tables simultaneously in memory (meaning, we
    // are only allocating tables for *at most* a [P C C E] but
    // it might be [P C E] or even [P E] depending on the
    // number of C sections in the structure (which depends on
    // the real unrolling amount).
    section_table_array.resize(numSectionStructures);
    unsigned sectionNo = 0;
    section_table_array[sectionNo++] = NULL;
    for (int p=0;p<numCoSectionStructures;p++) {
      section_table_array[sectionNo]) PartitionTables(Co);
      sectionNo++;
    }
    new (&section_table_array[sectionNo++])
      PartitionTables(E1);
    assert (sectionNo == section_structure_array.size());
#else
    // heterogeneous section_table_array; see LongTable case above 
    section_table_array.resize(numSectionStructures);
    for (unsigned p=0; p < section_table_array.size(); ++p) {
      section_table_array[p] = NULL;
    }
#endif
  } else if (tableOption == ZeroTable) {
    section_table_array.clear();
  } else {
    // assume NoTouchTable, so leave it alone.
    assert(tableOption == NoTouchTable);
  }
  
 #if 0 
  // FIXME  this is just commented out temporarily as refactoring moves code


 if (viterbiScore == true) {
    // Then we need to allocate the storage to keep the Viterbi options.
    // For now, it is unfortunately as long as the segment iself, even in the
    // island algorithm case, although the values area packed.

    // Once we include N-best list, then this is how it will effect these
    // data structures.
    const unsigned N_best = 1;

    // First do P
    P_partition_values.resize(N_best*section_structure_array[0].packer.packedLen());

    // Next C
    C_partition_values.resize(N_best
			      *section_structure_array[1].packer.packedLen()
			      *  (binaryViterbiFile ? 1 : numCoSectionTables) );

    // Next E
    E_partition_values.resize(N_best
			      *section_structure_array[numSectionStructures-1]
			      .packer.packedLen());

  }
#endif


  if (totalNumberSections != NULL)
    *totalNumberSections = (modifiedTemplateMaxUnrollAmount+3);

  return numUsableFrames;
}






// create the three junction trees for the basic sections.
void 
SectionScheduler::createSectionJunctionTrees(const string pStr) {
  infoMsg(IM::Giga,"Creating of P section JT\n");
  createSectionJunctionTree(gm_template.P,pStr);
  infoMsg(IM::Giga,"Creating of C section JT\n");
  createSectionJunctionTree(gm_template.C,pStr);
  infoMsg(IM::Giga,"Creating of E section JT\n");
  createSectionJunctionTree(gm_template.E,pStr);
  infoMsg(IM::Giga,"Done creating P,C,E section JTs\n");
}

// create the three junction forests for the basic sections
// from edges specified in trifile
void 
SectionScheduler::createSectionJunctionForests() {
  infoMsg(IM::Giga,"Creating of P section JT\n");
  createSectionJunctionForest(gm_template.P, string("default:P"));
  infoMsg(IM::Giga,"Creating of C section JT\n");
  createSectionJunctionForest(gm_template.C, string("default:C"));
  infoMsg(IM::Giga,"Creating of E section JT\n");
#if 1
  createSectionJunctionForest(gm_template.E, string("default:E"));
#else
  createSectionJunctionTree(gm_template.E,pStr);
#endif
  infoMsg(IM::Giga,"Done creating P,C,E section JTs\n");
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Support for building a tree from clique graph
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////






/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createSectionJunctionForest()
 *   Create a mini-junction forest from the cliques in the given section.
 *   This uses the edges specified by the message order given in the trifile
 *
 *   TODO: move this routine to a MaxClique class at some point.
 *   TODO: do this during triangulation time.
 *   TODO: this is section inference algorithm-specific
 *
 * Preconditions:
 *   The section must be instantiated with cliques 
 *
 * Postconditions:
 *   The cliques in the section are now such that they
 *   form a junction tree over cliques within that section.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   section.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void printUintSet(set<unsigned> s) {
  for (set<unsigned>::iterator it=s.begin(); it != s.end(); ++it) 
    printf(" %u", *it);
  printf("\n");
}
void 
SectionScheduler::createSectionJunctionForest(Section& section, string const &ia_name_and_section) {
  const unsigned num_max_cliques = section.cliques.size();

  infoMsg(IM::Giga,"Starting create JT\n");
  vector<set<unsigned> > clique_neighbors(num_max_cliques); // set of ith clique's neighbors
  vector< pair<unsigned, unsigned> > &msg_order = section.ia_message_order[ia_name_and_section];

  if (num_max_cliques == 0) {
    // Nothing to do.
    // This could happen if the section is empty which might occur
    // for empty P's and E's. C should never be empty.
    return;
  } else if (num_max_cliques == 1) {
    // then nothing to do
    infoMsg(IM::Giga,"Section has only one clique\n");
    return;
  } else {
    
    infoMsg(IM::Giga,"Section has %d cliques\n",num_max_cliques);
    if (section.ia_message_order.find(ia_name_and_section) == section.ia_message_order.end()) {
      error("ERROR: no message order found for inference architecture '%s'\n", 
	    ia_name_and_section.c_str());
    }
    for (unsigned msg_num=0; msg_num < msg_order.size(); ++msg_num) {
      unsigned i = msg_order[msg_num].first;
      unsigned j = msg_order[msg_num].second;
      clique_neighbors[i].insert(j);
      clique_neighbors[j].insert(i);
    }
  }

  // copy neighbor sets to vectors
  for (unsigned i=0; i < section.cliques.size(); ++i) {
    section.cliques[i].neighbors.clear();
    for (set<unsigned>::iterator it = clique_neighbors[i].begin();
	 it != clique_neighbors[i].end();
	 ++it)
    {
      section.cliques[i].neighbors.push_back(*it);
    }
  }
}









/*-
 *-----------------------------------------------------------------------
 * JunctionTree::createSectionJunctionTree()
 *   Create a mini-junction tree from the cliques in the given section.
 *   This uses Kruskal's greedy (but optimal) algorithm for MST generation.
 *
 *   TODO: move this routine to a MaxClique class at some point.
 *   TODO: do this during triangulation time.
 *
 * Preconditions:
 *   The section must be instantiated with cliques 
 *
 * Postconditions:
 *   The cliques in the section are now such that they
 *   form a junction tree over cliques within that section.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   section. Also initializes the Section's connected_components.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::createSectionJunctionTree(Section& section, const string junctionTreeMSTpriorityStr) {
  // TODO: junctionTreeMSTpriorityStr is shadowing the SectionScheduler::junctionTreeMSTpriorityStr
  //       member here. In SectionScheduler.h, the parameter is called pStr
  const unsigned numMaxCliques = section.cliques.size();

  infoMsg(IM::Giga,"Starting create JT\n");

  section.connected_components.clear();

  if (numMaxCliques == 0) {
    // Nothing to do.
    // This could happen if the section is empty which might occur
    // for empty P's and E's. C should never be empty.
    return;
  } else if (numMaxCliques == 1) {
    // then nothing to do
    infoMsg(IM::Giga,"Section has only one clique\n");
    section.connected_components.resize(1);
    section.connected_components[0].insert(0); // the single clique is the only component
    return;
  }

  infoMsg(IM::Giga,"Section has %d cliques\n",numMaxCliques);

  // Build the edge weight matrix representation of the clique graph.
  vector< Edge > edges;
  edges.reserve((numMaxCliques*(numMaxCliques-1))/2);

  for (unsigned i=0;i<numMaxCliques;i++) {
    for (unsigned j=i+1;j<numMaxCliques;j++) {
      set<RV*> sep_set;
      set_intersection(section.cliques[i].nodes.begin(),
		       section.cliques[i].nodes.end(),
		       section.cliques[j].nodes.begin(),
		       section.cliques[j].nodes.end(),
		       inserter(sep_set,sep_set.end()));
      Edge e;
      // define the edge
      e.clique1 = i; e.clique2 = j; 
      // !!!************************!!!
      // !!!** MUST DO THIS FIRST **!!!
      // !!!************************!!!
      // First push sep set size. To get a JT, we *MUST*
      // always choose from among the cliques that
      // have the largest intersection size.
      e.weights.push_back((double)sep_set.size());
      // now that the size is there, we have many other options.

      // All gets sorted in decreasing order, so that larger values
      // have priority. Thus, the remaining items we push back in
      // the case of ties.  
      //
      // *** Larger numbers are prefered. ***
      //
      // Options to include by priority of junctionTreeMSTpriorityStr:
      //   * D: number of deterministic nodes in separator 
      //   * E: number of deterministic nodes in union of cliques.
      //   * S: neg. weight of separator
      //   * U: neg. weight of union of cliques
      //   * V: neg. frame number variance in separator
      //   * W: neg. frame number variance in union
      //   * H: number of hidden nodes in separator
      //   * O: number of observed nodes in separator
      //   * L: number of hidden nodes in union
      //   * Q: number of observed nodes in union
      //
      //  Any can be preceeded by a '-' sign to flip effect. E.g., D-E-SU
      //
      // Default case: DSU
      //
      // TODO: Other ideas for this:
      //        a) maximize number of variables in same frame (or near each other) (like variance)
      //        b) minimize number of neighbors in each clique (i.e., 
      //           if cliques already have neighbors, choose the ones with fewer.
      //        c) integrate with RV value assignment to minimize
      //           the number of unassigned clique nodes (since they're
      //           iterated over w/o knowledge of any parents. If this
      //           ends up being a search, make this be offline, in with gmtkTriangulate
      // 

      float mult = 1.0;
      for (unsigned charNo=0;charNo< junctionTreeMSTpriorityStr.size(); charNo++) {
	const char curCase = toupper(junctionTreeMSTpriorityStr[charNo]);
	if (curCase == '-') {
	  mult = -1.0;
	  continue;
	}

	if (curCase == 'D') {

	  // push back number of deterministic nodes in
	  // the separator
	  set<RV*>::iterator it;
	  set<RV*>::iterator it_end = sep_set.end();
	  unsigned numDeterministicNodes = 0;
	  for (it = sep_set.begin(); it != it_end; it++) {
	    RV* rv = (*it);
	    if (rv->discrete() && RV2DRV(rv)->deterministic())
	      numDeterministicNodes++;
	  }
	  e.weights.push_back(mult*(double)numDeterministicNodes);

	} else if (curCase == 'E' || curCase == 'L' || curCase == 'Q' || curCase == 'W') {
	  // push back negative weight of two cliques together.
	  set<RV*> clique_union;
	  set_union(section.cliques[i].nodes.begin(),
		    section.cliques[i].nodes.end(),
		    section.cliques[j].nodes.begin(),
		    section.cliques[j].nodes.end(),
		    inserter(clique_union,clique_union.end()));

	  if (curCase == 'E') {
	    // push back number of deterministic nodes in
	    // the union
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = clique_union.end();
	    unsigned numDeterministicNodes = 0;
	    for (it = clique_union.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->discrete() && RV2DRV(rv)->deterministic())
		numDeterministicNodes++;
	    }
	    e.weights.push_back(mult*(double)numDeterministicNodes);
	  } else if (curCase == 'L') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = clique_union.end();
	    unsigned numHidden = 0;
	    for (it = clique_union.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->hidden())
		numHidden++;
	    }
	    e.weights.push_back(mult*(double)numHidden);
	  } else if (curCase == 'Q') {
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = clique_union.end();
	    unsigned numObserved = 0;
	    for (it = clique_union.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      if (rv->observed())
		numObserved++;
	    }
	    e.weights.push_back(mult*(double)numObserved);
	  } else if (curCase == 'W') {
	    // compute frame number variance in union, push back
	    // negative to prefer smalller frame variance (i.e.,
	    // connect things that on average are close to each
	    // other in time).
	    set<RV*>::iterator it;
	    set<RV*>::iterator it_end = clique_union.end();
	    double sum = 0;
	    double sumSq = 0;
	    for (it = clique_union.begin(); it != it_end; it++) {
	      RV* rv = (*it);
	      sum += rv->frame();
	      sumSq += rv->frame()*rv->frame();
	    }
	    double invsize = 1.0/(double)clique_union.size();
	    double variance = invsize*(sumSq - sum*sum*invsize);
	    e.weights.push_back(mult*(double)-variance);
	  }
	} else if (curCase == 'S') {

	  // push back negative weight of separator, to prefer
	  // least negative (smallest)  weight, since larger numbers
	  // are prefered.
	  e.weights.push_back(-(double)MaxClique::computeWeight(sep_set));

	  // printf("weight of clique %d = %f, %d = %f\n",
	  // i,section.cliques[i].weight(),
	  // j,section.cliques[j].weight());

	} else if (curCase == 'U') {

	  // push back negative weight of two cliques together.
	  set<RV*> clique_union;
	  set_union(section.cliques[i].nodes.begin(),
		    section.cliques[i].nodes.end(),
		    section.cliques[j].nodes.begin(),
		    section.cliques[j].nodes.end(),
		    inserter(clique_union,clique_union.end()));
	  e.weights.push_back(-(double)MaxClique::computeWeight(clique_union));
	} else if (curCase == 'V') {
	  // compute frame number variance in separator, push back
	  // negative to prefer smalller frame variance (i.e.,
	  // connect things that on average are close to each
	  // other in time).
	  set<RV*>::iterator it;
	  set<RV*>::iterator it_end = sep_set.end();
	  double sum = 0;
	  double sumSq = 0;
	  for (it = sep_set.begin(); it != it_end; it++) {
	    RV* rv = (*it);
	    sum += rv->frame();
	    sumSq += rv->frame()*rv->frame();
	  }
	  if (sep_set.size() == 0) {
	    e.weights.push_back(mult*(double)-FLT_MAX);
	  } else {
	    double invsize = 1.0/(double)sep_set.size();
	    double variance = invsize*(sumSq - sum*sum*invsize);
	    e.weights.push_back(mult*(double)-variance);
	  }
	} else if (curCase == 'H') {
	  set<RV*>::iterator it;
	  set<RV*>::iterator it_end = sep_set.end();
	  unsigned numHidden = 0;
	  for (it = sep_set.begin(); it != it_end; it++) {
	    RV* rv = (*it);
	    if (rv->hidden())
	      numHidden++;
	  }
	  e.weights.push_back(mult*(double)numHidden);
	} else if (curCase == 'O') {
	  set<RV*>::iterator it;
	  set<RV*>::iterator it_end = sep_set.end();
	  unsigned numObserved = 0;
	  for (it = sep_set.begin(); it != it_end; it++) {
	    RV* rv = (*it);
	    if (rv->observed())
	      numObserved++;
	  }
	  e.weights.push_back(mult*(double)numObserved);
	} else {
	  error("ERROR: Unrecognized junction tree clique sort order letter '%c' in string '%s'\n",
		curCase, junctionTreeMSTpriorityStr.c_str());
	}
	mult = 1.0;
      }

      // add the edge.
      edges.push_back(e);
      if (IM::messageGlb(IM::Giga)) {
	infoMsg(IM::Giga,"Edge (%d,%d) has sep size %.0f, ",
		i,j,
		e.weights[0]);
	for (unsigned charNo=0;charNo< junctionTreeMSTpriorityStr.size(); charNo++) {
	  const char curCase = toupper(junctionTreeMSTpriorityStr[charNo]);
	  infoMsg(IM::Giga,"%c,weight[%d] = %f, ",curCase,charNo+1,e.weights[charNo+1]);
	}
	infoMsg(IM::Giga,"\n");
      }
    }
  }
  // sort in decreasing order by edge weight which in this
  // case is the sep-set size.
  sort(edges.begin(),edges.end(),EdgeCompare());

  // Run max spanning tree to construct JT from set of cliques.
  // This is basically Kruskal's algorithm, but without using the
  // fast data structures (it doesn't need to be that fast
  // since it is run one time per section, for all inference
  // runs of any length.
  
  // Create a vector of sets, corresponding to the trees associated
  // with each clique. Non-empty set intersection corresponds to
  // tree overlap. Each set contains the clique indices in teh set.
  vector < set<unsigned> >  findSet;
  findSet.resize(numMaxCliques);
  for (unsigned i=0;i<numMaxCliques;i++) {
    set<unsigned> iset;
    iset.insert(i);
    findSet[i] = iset;
  }
  
  unsigned joinsPlusOne = 1;
  for (unsigned i=0;i<edges.size();i++) {
    infoMsg(IM::Giga,"Edge %u has sep size %.0f\n",
	    i,
	    edges[i].weights[0]);
    
    if (edges[i].weights[0] == 0.0) {
      if (IM::messageGlb(IM::High)) {
	// there is no way to know the difference here if the
	// graph is non-triangualted or is simply disconnected
	// (which is ok). A non-triangulated graph might have
	// resulted from the user editing the trifile, but we
	// presume that MCS has already checked for this when
	// reading in the trifiles. We just issue an informative
	// message just in case.
	printf("NOTE: junction tree creation not joining two cliques (%d and %d) with size 0 set intersection. Either disconnected (which is ok) or non-triangulated (which is bad) graph.\n",
	       edges[i].clique1,edges[i].clique2);
	// TODO: print out two cliques that are trying to be joined.
	printf("Clique %d: ",edges[i].clique1);
	section.cliques[edges[i].clique1].printCliqueNodes(stdout);
	printf("Clique %d: ",edges[i].clique2);
	section.cliques[edges[i].clique2].printCliqueNodes(stdout);
      }
    } else {
      // Only use non-empty sep set edges to construct junction forest
      // to avoid connecting a disconnected graph 

      set<unsigned>& iset1 = findSet[edges[i].clique1];
      set<unsigned>& iset2 = findSet[edges[i].clique2];
      if (iset1 != iset2) {
	// merge the two sets
	set<unsigned> new_set;
	set_union(iset1.begin(),iset1.end(),
		  iset2.begin(),iset2.end(),
		  inserter(new_set,new_set.end()));
	// make sure that all members of the set point to the
	// new set.
	set<unsigned>::iterator ns_iter;
	for (ns_iter = new_set.begin(); ns_iter != new_set.end(); ns_iter ++) {
	  const unsigned clique = *ns_iter;
	  findSet[clique] = new_set;
	}

        infoMsg(IM::Giga,"Joining cliques %d and %d (edge %d) with intersection size %.0f\n",
		edges[i].clique1,edges[i].clique2,i,edges[i].weights[0]);
        section.cliques[edges[i].clique1].neighbors.push_back(edges[i].clique2);
        section.cliques[edges[i].clique2].neighbors.push_back(edges[i].clique1);
      }
      // Early termination optimization: there are (n)(n-1)/2 edges to iterate 
      //   over, but a spanning forest can use at most n-1.  I think it's safe to
      //   count zero weight edges here - actually, if the edges are sorted with
      //   sep set size as the major key, I think we could break on 
      //   weights[0] == 0.0  -RR
      if (++joinsPlusOne == numMaxCliques)
	break;
    }
  }
  
  // Initialize connected components using left-over Kruskal union-find data

  // For clique # i, findSet[i] is the set of cliques in the same
  // connected component as clique # i. If clique #s i and j are in
  // the same connected component, findSet[i] == findSet[j], so we
  // want to find the set of unique findSet values.

  // visited[i] = if we already know clique #i's connected component
  vector<bool> visited(numMaxCliques, false); 
  unsigned num_components = 0;
  for (unsigned i = 0; i < numMaxCliques; ++i) {
    set<unsigned> &i_set = findSet[i];
    set<unsigned>::iterator it = i_set.begin(); // first clique # in this connected component
    if (! visited[*it] ) {
      // we've discovered a new component
      section.connected_components.push_back(set<unsigned>()); // grow the connected_components vector
      for ( ; it != i_set.end(); ++it) {
	section.connected_components[num_components].insert(*it);
	visited[*it] = true;
      }
      num_components += 1;
    }
  }

  printf("connected components:\n");
  for (unsigned i=0; i < num_components; ++i) {
    printf("  %u: ", i);
    set<unsigned> &cc = section.connected_components[i];
    for (set<unsigned>::iterator it = cc.begin(); it != cc.end(); ++it) {
      printf("%u {", *it);
      printRVSet(stdout, section.cliques[*it].nodes, false);
      printf("}  ");
    }
    printf("\n");
  }

}



/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::create_base_sections()
 *   sets up the internal data structures and creates a junction tree corresponding
 *   to the graph with C' unrolled k=1 times, (i.e., we have
 *   a graph that has (1) copies of C'). 
 *
 *   This routine should only be called once to set up the base
 *   sections (P',C',E'), from which all future unrolls will take
 *   place, modifications, and instantiations take place.
 *
 *   Sets things up so that inference can take place.
 *
 *   Unrolling is in terms of modifid (i.e., C') sections (so that
 *   if S=1, k corresponds to frames). See GMTemplate::computeUnrollParameters
 *   for documentation on how this relates to frames, and the basic and modified
 *   template.
 *
 * Preconditions:
 *   The sections must be validly instantiated with cliques 
 *
 * Postconditions:
 *   The sections are such that they now know who they are connected to on
 *   the right and left.
 *
 * Side Effects:
 *   Modifies all neighbors variables within the cliques within the
 *   section.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::create_base_sections()
{

  // first create the unrolled set of random variables corresponding
  // to this JT. Unroll corresponding to 1 [P' C' E']. (see GMTemplate::computeUnrollParameters
  // for more documentation on this).

  unsigned M = gm_template.maxNumChunksInBoundary();
  unsigned S = gm_template.chunkSkip();

  fp.unroll(M + S - 1,
	    section_unrolled_rvs,section_ppf);

  vector< set <RV*> > empty;

  // copy P section 
  new (&P1) JT_Partition(gm_template.P, 0*S*fp.numFramesInC(),
			 empty, 0*S*fp.numFramesInC(),
			 gm_template.PCInterface_in_P, 0*S*fp.numFramesInC(),
			 section_unrolled_rvs,section_ppf);

  P1_message_order = P1.ia_message_order[string("default:P")];
  // copy E section
assert(gm_template.E.cliques.size() > 0);
  new (&E1) JT_Partition(gm_template.E, 0*S*fp.numFramesInC(),
			 gm_template.CEInterface_in_E, 0*S*fp.numFramesInC(),
			 empty, 0*S*fp.numFramesInC(),
			 section_unrolled_rvs,section_ppf);
  E1_message_order = E1.ia_message_order[string("default:E")];

  if (gm_template.usesLeftInterface()) {
    // left interface case
    // Template has a (P)(C)(CE) = (P')(C')(E')

    // An E' in the LI case must be E' = [ C E ], so there is never an
    // empty E here (since C is never empty).
    assert ( E1.cliques.size() > 0 );

    // there might be an empty P1, though, so we need to check for that.
    if (P1.cliques.size() > 0) {
      // neither P1 nor E1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    } else {
      // P1 is empty. For Co's left interface, we use its right
      // interface since there is no interface to P1. The reason why
      // this is valid is that E' = [C E], and E can't connect to C1
      // in [C1 C2 E], so C2's li to C1 is the same as E's li to C' in
      // [C' E']. We do need to shift E'ls li to C' to the left by S
      // though.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.CEInterface_in_C,-1*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    }
  } else {
    // right interface case, symmetric to the left interface case above.
    // Right interface case.
    // Template has a (PC)(C)(E) = (P')(C')(E')

    // An P' in the RI case must be P' = [ P C ], so
    // there is never an empty P' here.
    assert ( P1.cliques.size() > 0 );

    // there might be an empty E1, though, so we need to check for that.

    if (E1.cliques.size() > 0) {
      // neither E1 nor P1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    } else {
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.PCInterface_in_C,1*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    }
  }
  Co_message_order = Co.ia_message_order[string("default:C")];
}

// This version supports mean field approximation inference architecture
void
SectionScheduler::create_base_sections_mfa()
{

  // first create the unrolled set of random variables corresponding
  // to this JT. Unroll corresponding to 1 [P' C' E']. (see GMTemplate::computeUnrollParameters
  // for more documentation on this).

  unsigned M = gm_template.maxNumChunksInBoundary();
  unsigned S = gm_template.chunkSkip();

  fp.unroll(M + S - 1,
	    section_unrolled_rvs,section_ppf);

  vector< set <RV*> > empty;

  // copy P section 
  new (&P1) JT_Partition(gm_template.P, 0*S*fp.numFramesInC(),
			 empty, 0*S*fp.numFramesInC(),
			 gm_template.PCInterface_in_P, 0*S*fp.numFramesInC(),
			 section_unrolled_rvs,section_ppf);

  // copy E section
  new (&E1) JT_Partition(gm_template.E, 0*S*fp.numFramesInC(),
			 gm_template.CEInterface_in_E, 0*S*fp.numFramesInC(),
			 empty, 0*S*fp.numFramesInC(),
			 section_unrolled_rvs,section_ppf);

  if (gm_template.usesLeftInterface()) {
    // left interface case
    // Template has a (P)(C)(CE) = (P')(C')(E')

    // An E' in the LI case must be E' = [ C E ], so there is never an
    // empty E here (since C is never empty).
    assert ( E1.cliques.size() > 0 );

    // there might be an empty P1, though, so we need to check for that.
    if (P1.cliques.size() > 0) {
      // neither P1 nor E1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    } else {
      // P1 is empty. For Co's left interface, we use its right
      // interface since there is no interface to P1. The reason why
      // this is valid is that E' = [C E], and E can't connect to C1
      // in [C1 C2 E], so C2's li to C1 is the same as E's li to C' in
      // [C' E']. We do need to shift E'ls li to C' to the left by S
      // though.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.CEInterface_in_C,-1*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    }
  } else {
    // right interface case, symmetric to the left interface case above.
    // Right interface case.
    // Template has a (PC)(C)(E) = (P')(C')(E')

    // An P' in the RI case must be P' = [ P C ], so
    // there is never an empty P' here.
    assert ( P1.cliques.size() > 0 );

    // there might be an empty E1, though, so we need to check for that.

    if (E1.cliques.size() > 0) {
      // neither E1 nor P1 are empty.
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.CEInterface_in_C,0*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    } else {
      new (&Co) JT_Partition(gm_template.C,0*S*fp.numFramesInC(),
			     gm_template.PCInterface_in_C,0*S*fp.numFramesInC(),			   
			     gm_template.PCInterface_in_C,1*S*fp.numFramesInC(),
			     section_unrolled_rvs,section_ppf);
    }
  }
}

// routine to find the interface cliques of the sections
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::computeSectionInterfaces()
 *   This version finds the vairous "interface" between the sections
 *   in a graph, but does so for all possible unrollings by doing
 *   one unrolling. It fills in the P_ri_to_C, C_li_to_P variables
 *   and is thus necessary to be called before any true junction tree
 *   creation.
 *
 *
 * Preconditions:
 *   The gm_template sections must be validly instantiated with cliques.
 *
 * Postconditions:
 *   The member variables P_ri_to_C etc. are filled in.
 *
 * Side Effects:
 *   chnages group of member variables P_ri_to_C, C_li_to_P, etc.
 *
 * Results:
 *   Modifies member variables as above.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::computeSectionInterfaces() {

  // set up the base sections
  create_base_sections();


#if 1
  P_ri_to_C = P1.section_ri;
  P_ri_to_C_size = P_ri_to_C.size(); // gm_template.PCInterface_in_P.size();
  E_li_to_C = E1.section_li;
  
  C_li_to_P = Co.section_li;
  C_li_to_C = C_li_to_P;
  C_ri_to_E = Co.section_ri;
  C_ri_to_C = C_ri_to_E;

  // sanity check
  assert (C_ri_to_C == C_ri_to_E);
  
  C_ri_to_C_size = C_ri_to_C.size(); //gm_template.CEInterface_in_C.size();

  
  set<unsigned> potential_E_roots;
  for (unsigned i=0; i < E1.cliques.size(); ++i) {
    potential_E_roots.insert(i);
  }
  for (vector<pair<unsigned,unsigned> >::iterator it = E1_message_order.begin();
       it != E1_message_order.end();
       ++it)
  {
    potential_E_roots.erase((*it).first); // roots cannot be message source
  }
  for (set<unsigned>::iterator it = potential_E_roots.begin();
       it != potential_E_roots.end();
       ++it)
  {
    E_root_clique.push_back(*it);
  }

#else
  // Use base sections to find the various interface cliques.
  bool P_riCliqueSameAsInterface;
  bool E_liCliqueSameAsInterface;
  P1.findRInterfaceClique(P_ri_to_C,P_riCliqueSameAsInterface,interfaceCliquePriorityStr);
  P_ri_to_C_size = gm_template.PCInterface_in_P.size();
  E1.findLInterfaceClique(E_li_to_C,E_liCliqueSameAsInterface,interfaceCliquePriorityStr);

  // Note that here, we do the same for both left and right interface,
  // but we break it out into two cases for clarity.
  if (gm_template.usesLeftInterface()) {
    // left interface case

    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Co.findLInterfaceClique(C_li_to_C,Co_liCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_li_to_P = C_li_to_C;
    Co.findRInterfaceClique(C_ri_to_E,Co_riCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_ri_to_C = C_ri_to_E;

    // sanity check
    assert (C_ri_to_C == C_ri_to_E);

    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Co_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Co_riCliqueSameAsInterface && E_liCliqueSameAsInterface;

  } else { // right interface

    bool Co_liCliqueSameAsInterface;
    bool Co_riCliqueSameAsInterface;

    Co.findLInterfaceClique(C_li_to_P,Co_liCliqueSameAsInterface,interfaceCliquePriorityStr);
    Co.findRInterfaceClique(C_ri_to_C,Co_riCliqueSameAsInterface,interfaceCliquePriorityStr);
    C_li_to_C = C_li_to_P;
    C_ri_to_E = C_ri_to_C;

    // sanity check
    assert (C_li_to_C == C_li_to_P);

    P_to_C_icliques_same =
      P_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;

    C_to_C_icliques_same =
      Co_riCliqueSameAsInterface && Co_liCliqueSameAsInterface;    

    C_to_E_icliques_same =
      Co_riCliqueSameAsInterface && E_liCliqueSameAsInterface;
  }
  C_ri_to_C_size = gm_template.CEInterface_in_C.size();

  // TODO: see if it is possible to choose a better root for E.  
  // TODO: make command line heuristics for choosing E-root-clique as well.
  // E order, clique 0 could be choosen as root arbitrarily.  
  // E_root_clique = 0;
  // E_root_clique = E1.cliqueWithMinWeight();
  E_root_clique.resize(1);
  E_root_clique[0] = E1.cliqueWithMaxWeight();
  // If this is updated, need also to update in all other places
  // the code, search for string "update E_root_clique"
#endif
}

// This version supports mean field approximation inference architecture
void 
SectionScheduler::computeSectionInterfacesMFA() {


  // FIXME - put these somewhere
#if 0
  P1_message_order = P1.ia_message_order[string("default:P")];
  E1_message_order = E1.ia_message_order[string("default:E")];
  Co_message_order = Co.ia_message_order[string("default:C")];
#endif

  // set up the base sections
  create_base_sections_mfa();

  printf("Finding interface clique(s) amoung:\n");
  bool P_riCliqueSameAsInterface;
  bool E_liCliqueSameAsInterface;
printf("P -> C interface:\n");
  P1.findRInterfaceClique(P_ri_to_C,P_riCliqueSameAsInterface,interfaceCliquePriorityStr);
printf("C <- E interface:\n");
  E1.findLInterfaceClique(E_li_to_C,E_liCliqueSameAsInterface,interfaceCliquePriorityStr);


  twiddle according to use left/right interface. see JunctionTree::computePartitionInterface()
  pick roots for P, C, E sub JTs - not necessarily interface cliques as there might be no temporal edges
		      - here, or in computeMessageOrders ?

  P_ri_to_C = P1.section_ri;
  P_ri_to_C_size = P_ri_to_C.size(); // gm_template.PCInterface_in_P.size();
  E_li_to_C = E1.section_li;
  
  C_li_to_P = Co.section_li;
  C_li_to_C = C_li_to_P;
  C_ri_to_E = Co.section_ri;
  C_ri_to_C = C_ri_to_E;


  // FIXME - need to pick an E root for each connected component, no right interface
  // sanity check
  assert (C_ri_to_C == C_ri_to_E);
  
  C_ri_to_C_size = C_ri_to_C.size(); //gm_template.CEInterface_in_C.size();

  
  set<unsigned> potential_E_roots;
  for (unsigned i=0; i < E1.cliques.size(); ++i) {
    potential_E_roots.insert(i);
  }
  for (vector<pair<unsigned,unsigned> >::iterator it = E1_message_order.begin();
       it != E1_message_order.end();
       ++it)
  {
    potential_E_roots.erase((*it).first); // roots cannot be message source
  }
  for (set<unsigned>::iterator it = potential_E_roots.begin();
       it != potential_E_roots.end();
       ++it)
  {
    E_root_clique.push_back(*it);
  }
}




/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createFactorCliques()
 *   place factorClique in either P1, Co, or E whichever	 
 *   is the left most that contains all the variables.	 
 *   this is thus optimized for left-to-right inference.
 *
 *
 * Preconditions:
 *   Should be called only from createFactorClique.
 *
 * Postconditions:
 *   The factor cliques have been inserted in the appropraite JT_Partitions
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::insertFactorClique(FactorClique& factorClique,FactorInfo& factor)
{
  count_iterator<set <RV*> > res;
  // first try P1
  set_intersection(factorClique.nodes.begin(),factorClique.nodes.end(),
		   P1.nodes.begin(),P1.nodes.end(),
		   res);
  if (res.count() == factorClique.nodes.size()) {
    // then fully contained in P1
    infoMsg(IM::Giga,"insertFactorClique: inserting factor %s(%d) into section %s\n",
	    factor.name.c_str(),factor.frame,P1_n);
    P1.factorCliques.push_back(factorClique);
  } else {
    // try Co
    set_intersection(factorClique.nodes.begin(),factorClique.nodes.end(),
		     Co.nodes.begin(),Co.nodes.end(),
		     res);
    if (res.count() == factorClique.nodes.size()) {
      // then fully contained in P1
//Co ?
      infoMsg(IM::Giga,"insertFactorClique: inserting factor %s(%d) into section %s\n",
	      factor.name.c_str(),factor.frame,Co_n);
      Co.factorCliques.push_back(factorClique);
    } else {
      // try E1
      set_intersection(factorClique.nodes.begin(),factorClique.nodes.end(),
		       E1.nodes.begin(),E1.nodes.end(),
		       res);
      if (res.count() == factorClique.nodes.size()) {
	// then fully contained in P1
//E1 ?
	infoMsg(IM::Giga,"insertFactorClique: inserting factor %s(%d) into section %s\n",
		factor.name.c_str(),factor.frame,E1_n);
	E1.factorCliques.push_back(factorClique);
      } else {
	error("INTERNAL ERROR: factor %s(%d) defined at %s:%d contains nodes that does not live in any section\n",
	      factor.name.c_str(),factor.frame,
	      factor.fileName.c_str(),factor.fileLineNumber);
	assert(0); // this should never happen.
      }
    }
  }
}





// routine to create the factors in the appropriate sections
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createFactorCliques()
 *   Create the factor cliques (i.e., hard/soft) constraints to
 *   be part of the sections and maxcliques.
 * 
 *
 * Preconditions:
 *   Sections must be instantiated, and interface cliques
 *   must have been computed (i.e., computeSectionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   The factor cliques have been created in the three main JT_Partitions
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::createFactorCliques() {
  unsigned M = gm_template.maxNumChunksInBoundary();
  unsigned S = gm_template.chunkSkip();

  // do this for all original C sections.
  // since this is the basic section, we have
  // M+S copies of the chunk.
  const unsigned numChunkCopies = M + S;

  for (unsigned factorNo=0;factorNo < fp.factorList.size(); factorNo++) {
    // find a home for the current factor, in either
    // P1, Co, or E1, being mindful of the partial
    // unrolling that potentially has already taken
    // place in producing P1, Co, and E1.

    FactorInfo& factor = fp.factorList[factorNo];

    if (fp.frameInTemplateP(factor.frame)) {
      const unsigned offset = 0;

      FactorClique factorClique(factor,
				section_unrolled_rvs,section_ppf,
				offset);

      insertFactorClique(factorClique,factor);

    } else if (fp.frameInTemplateC(factor.frame)) {

      for (unsigned i=0;i<numChunkCopies;i++) {
	const unsigned offset = 
	  i*fp.numFramesInC();

	FactorClique factorClique(factor,
				  section_unrolled_rvs,section_ppf,
				  offset);
	insertFactorClique(factorClique,factor);
      }
    } else {
      assert (fp.frameInTemplateE(factor.frame));
      const unsigned offset = 
	(numChunkCopies-1)*fp.numFramesInC();

      FactorClique factorClique(factor,
				section_unrolled_rvs,section_ppf,
				offset);
      insertFactorClique(factorClique,factor);
    }
  }
}


// routine to find the interface cliques of a section
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::computeSectionInterface()
 *   Finds the "interface" between the two sections. The interface are
 *   the two cliques in each interface that have maximum intersection.
 *   These are the two cliques that are joined in a junction tree.
 *   It is assumed that section 1 is just to the "left" of section 2.
 *   This does an N*M algorithm, where N (resp. M) is the number of cliques 
 *   in one (resp. the other) section
 *
 *
 * Preconditions:
 *   The sections must be validly instantiated with cliques, and section 1
 *   must be just to the left of section 2 (i.e., valid combinations
 *   include (P,C), (C,E), or (C1,C2). 
 *
 * Postconditions:
 *   The sections are such that they now know who they are connected to on
 *   the right and left.
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   The left sections right interface clique and the right sections
 *   left interface clique are returned via the section1_ric and section2_lic variables.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::computeSectionInterface(JT_Partition& section1,
					// section 1's right interface clique
					  unsigned int& section1_ric,
					  JT_Partition& section2,
					// section 2's left interface clique
					  unsigned int& section2_lic,
					// are the two interface cliques identical
					  bool& icliques_same)
{
  const unsigned numMaxCliques1 = section1.cliques.size();
  const unsigned numMaxCliques2 = section2.cliques.size();

  // run a dumb n^2 algorithm to find the two cliques in each
  // section that have maximum overlap. Because the random variable
  // set is not guaranteed to be the same at this point, we need to
  // compute intersection by (name,frame) equivalance.

  infoMsg(IM::Giga,"Starting create section interface\n");

  unsigned max_sep_set_size = 0;
  unsigned section1_clique=0,section2_clique=0;
  unsigned section1_clique_size=0,section2_clique_size=0;
  for (unsigned i=0;i<numMaxCliques1;i++) {
    for (unsigned j=0;j<numMaxCliques2;j++) {

      set< RV*>::iterator it1;
      const set< RV*>::iterator it1_end = section1.cliques[i].nodes.end();
      set< RV*>::iterator it2;
      const set< RV*>::iterator it2_end = section2.cliques[j].nodes.end();

      unsigned sep_set_size = 0;
      for (it1 = section1.cliques[i].nodes.begin(); it1 != it1_end ; it1++) {
	for (it2 = section2.cliques[j].nodes.begin(); it2 != it2_end ; it2++) {
	  RV* rv1 = (*it1);
	  RV* rv2 = (*it2);
	  infoMsg(IM::Giga,"comparing %s(%d) with %s(%d)\n",
		  rv1->name().c_str(),rv1->frame(),		  
		  rv2->name().c_str(),rv2->frame());
	  if (rv1->frame() == rv2->frame() && rv1->name() == rv2->name()) {
	    sep_set_size++;
	  }
	}
      }
      infoMsg(IM::Giga,"clique %d of section1 and clique %d of section2 has overlap %d\n",i,j,sep_set_size);
      if (sep_set_size > max_sep_set_size) {
	section1_clique = i;
	section1_clique_size = section1.cliques[i].nodes.size();
	section2_clique = j;
	section2_clique_size = section2.cliques[j].nodes.size();
	max_sep_set_size = sep_set_size;
      }
    }
  }
  infoMsg(IM::Giga,"clique %d of section1 (size %d) and clique %d of section2 (size %d) has max overlap of %d\n",
	  section1_clique,section1_clique_size,
	  section2_clique,section2_clique_size,
	  max_sep_set_size);
  section1_ric = section1_clique;
  section2_lic = section2_clique;
  if ((section1_clique_size == section2_clique_size)
      &&
      (section1_clique_size == max_sep_set_size)) {
    // then the interface cliques on both sides of the
    // section are the same, so one can be removed.
    icliques_same = true;
  } else {
    icliques_same = false;
  }
}


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::prepareForNextInferenceRound()
 *   does a bit if setup for next inference round.
 *
 * Preconditions:
 *   The sections must be validly instantiated with cliques, and
 *   the routine assignRVsToCliques() must have been called.
 *
 * Postconditions:
 *   initial values are re-initialized.
 *
 * Side Effects:
 *   modifies a few variables in sections.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::prepareForNextInferenceRound()
{
  for (unsigned c=0;c<P1.cliques.size();c++)
    P1.cliques[c].prepareForNextInferenceRound();
  for (unsigned c=0;c<Co.cliques.size();c++)
    Co.cliques[c].prepareForNextInferenceRound();
  for (unsigned c=0;c<E1.cliques.size();c++)
    E1.cliques[c].prepareForNextInferenceRound();
}

// root the JT
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createDirectedGraphOfCliques()
 *   Turns undirected graph into DAG of cliques, for all sections.
 * 
 *
 * Preconditions:
 *   Sections must be instantiated, and interface cliques
 *   must have been computed (i.e., computeSectionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   The member variables P_to_C_message_order (etc.) are created
 *   with the order of message passing.
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
printChildren(char name, vector<MaxClique> &cliques) {
  printf("%c children:\n", name);
  for (unsigned i=0; i < cliques.size(); ++i) {
    printf("  %u:", i);
    for (unsigned j=0; j < cliques[i].children.size(); ++j)
      printf(" %u", cliques[i].children[j]);
    printf("\n");
  }
}

void 
SectionScheduler::createDirectedGraphOfCliques() {
  createDirectedGraphOfCliques(P1, P_ri_to_C);
  if (gm_template.usesLeftInterface()) {
    createDirectedGraphOfCliques(Co, C_ri_to_E);
  } else { // right interface
    createDirectedGraphOfCliques(Co, C_ri_to_C);
  }
  createDirectedGraphOfCliques(E1, E_root_clique);
}


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createDirectedGraphOfCliques()
 * SectionScheduler::createDirectedGraphOfCliquesRecurse()
 *   From the root that is given as an argument, turn the
 *   undirected graph of cliques JT into a directed graph
 *   of cliques rooted at root. Do this by assigning
 *   to each clique's children member variable.
 *
 * Preconditions:
 *   Sections must be instantiated, and interface cliques
 *   must have been computed (i.e., computeSectionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   children member variables in cliques are set.
 *
 * Side Effects:
 *    none
 *
 * Results:
 *    None, but passed back in arguments.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::createDirectedGraphOfCliques(JT_Partition& part,
					   const unsigned root)
{
  // check for empty partitions, possible for P or E.
  if (part.cliques.size() == 0)
    return;

  // make sure none have been so far visited
  vector< bool > visited(part.cliques.size());
  for (unsigned i=0;i<part.cliques.size();i++) {
    visited[i] = false;
    part.cliques[i].children.clear();
  }
  createDirectedGraphOfCliquesRecurse(part,root,visited);
}
void 
SectionScheduler::createDirectedGraphOfCliques(JT_Partition& section, vector<unsigned> const &roots) {
  // check for empty sections, possible for P or E.
  if (section.cliques.size() == 0)
    return;

  // make sure none have been so far visited
  vector< bool > visited(section.cliques.size());
  for (unsigned i=0;i<section.cliques.size();i++) {
    visited[i] = false;
    section.cliques[i].children.clear();
  }

  // The undirected graph is built from the message order specified in 
  // the trifile inference architecture, and thus may not be connected.
  // In GMTK 1.X, the undirected graph was built with the MST code 
  // in JunctionTree::createPartitionJunctionTrees(), which treats all
  // the cliques as a complete graph ensuring the undirected graph is
  // connected. But in GMTK 2.X, the undirected graph may not be connected, 
  // necessitating a loop over the roots to produce a directed junction forest.

  for (unsigned i=0; i < roots.size(); ++i) {
    if (!visited[roots[i]]) {
      createDirectedGraphOfCliquesRecurse(section,roots[i],visited);
    }
  }
  // Furthermore, there may be connected components in the undirected 
  // graph that do not contain any "root" (interface) cliques, so we 
  // still have to look for any unvisited cliques.

  for (unsigned i=0;i<section.cliques.size();i++) {
    if (!visited[i]) {
      createDirectedGraphOfCliquesRecurse(section,i,visited);
    }
  }


}

void
SectionScheduler::createDirectedGraphOfCliquesRecurse(JT_Partition& section,
						  const unsigned root,
						  vector< bool >& visited)
{
  visited[root] = true;
  for (unsigned childNo=0;
       childNo<section.cliques[root].neighbors.size();
       childNo++) {
    const unsigned child = section.cliques[root].neighbors[childNo];
    if (visited[child])
      continue;
    section.cliques[root].children.push_back(child);
    createDirectedGraphOfCliquesRecurse(section,child,visited);
  }
}


// Assign probability giving random variables to cliques (i.e.,
// these are assigned only to cliques such that the random variables
// and *all* their parents live in the clique, plus some other
// criterion in order to make message passing as efficient as
// possible).
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::assignRVsToCliques()
 *    Assign to the cliques the random variables that will
 *    be used to generate the tables for each clique. This is
 *    mainly a stub routine to operate on all sections. The
 *    real code is in the other routines below.
 *
 * Preconditions:
 *     createSectionJunctionTrees(), computeSectionInterfaces(), and
 *     createDirectedGraphOfCliques() must have been called before 
 *     this can be called.
 *
 * Postconditions:
 *     Each of the sections so specified in the code below
 *     are now such that their cliques have their 'assignedNodes' 
 *     member filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::
assignRVsToCliques(const char* varSectionAssignmentPrior, const char *varCliqueAssignmentPrior) {
#if 1
#else
  if (P1.nodes.size() > 0) {
    infoMsg(IM::Giga,"assigning rvs to P1 section\n");
    assignRVsToCliques(P1_n,P1,/*P_ri_to_C*/ P1.connected_jt_root,varSectionAssignmentPrior,varCliqueAssignmentPrior);
    // accumulate what occured in P1 so Co can use it. 
#if 1
    unionRVs(P1.cliques[P1.connected_jt_root].cumulativeAssignedNodes,
	     P1.cliques[P1.connected_jt_root].assignedNodes,
	     Co.cliques[C_li_to_P].cumulativeAssignedNodes);
    unionRVs(P1.cliques[P_ri_to_C].cumulativeAssignedProbNodes,
	     P1.cliques[P_ri_to_C].assignedProbNodes,
	     Co.cliques[C_li_to_P].cumulativeAssignedProbNodes);
#else
    unionRVs(P1.cliques[P_ri_to_C].cumulativeAssignedNodes,
	     P1.cliques[P_ri_to_C].assignedNodes,
	     Co.cliques[C_li_to_P].cumulativeAssignedNodes);
    unionRVs(P1.cliques[P_ri_to_C].cumulativeAssignedProbNodes,
	     P1.cliques[P_ri_to_C].assignedProbNodes,
	     Co.cliques[C_li_to_P].cumulativeAssignedProbNodes);
#endif
    // printf("P1's cum ass prob:");
    // printRVSet(stdout,P1.cliques[P_ri_to_C].cumulativeAssignedProbNodes);
    // printf("P1's ass prob:");
    // printRVSet(stdout,P1.cliques[P_ri_to_C].assignedProbNodes);
    // printf("Co's cum ass prob:");
    // printRVSet(stdout,Co.cliques[C_li_to_P].cumulativeAssignedProbNodes);

  }
  infoMsg(IM::Max,"assigning rvs to Co section\n");
  assignRVsToCliques(Co_n,Co,/*C_ri_to_C*/ Co.connected_jt_root,varSectionAssignmentPrior,varCliqueAssignmentPrior);

  if (E1.nodes.size() > 0) {
    infoMsg(IM::Max,"assigning rvs to E section\n");
    // accumulate what occured in P1,Co so E1 can use it. 
    unionRVs(Co.cliques[C_ri_to_E].cumulativeAssignedNodes,
	     Co.cliques[C_ri_to_E].assignedNodes,
	     E1.cliques[E_li_to_C].cumulativeAssignedNodes);
    unionRVs(Co.cliques[C_ri_to_E].cumulativeAssignedProbNodes,
	     Co.cliques[C_ri_to_E].assignedProbNodes,
	     E1.cliques[E_li_to_C].cumulativeAssignedProbNodes);
    assignRVsToCliques(E1_n,E1,/*E_root_clique*/ E1.connected_jt_root,varSectionAssignmentPrior,varCliqueAssignmentPrior);
  }

  // lastly, check to make sure all nodes have been assigned to to
  // give probability to one clique. This should be done in case the
  // user edited the trifile, which is one possible reason for this to
  // end in error.
  set <RV*> allNodes;
  unionRVs(P1.nodes,Co.nodes,allNodes);
  unionRVs(Co.nodes,E1.nodes,allNodes,true);
  set <RV*> allAssignedProbNodes;
  if (E1.cliques.size() > 0)
    unionRVs(E1.cliques[E_root_clique].cumulativeAssignedProbNodes,
	     E1.cliques[E_root_clique].assignedProbNodes,
	     allAssignedProbNodes);
  else 
    unionRVs(Co.cliques[C_ri_to_E].cumulativeAssignedProbNodes,
	     Co.cliques[C_ri_to_E].assignedProbNodes,
	     allAssignedProbNodes);

  set <RV*> nodesThatGiveNoProb;
// 153: OK nodesThatGiveNoProb printed in error message
// 153:    though we could use count_iterator to be efficient in
// 153:    the common case, and compute the intersection in the
// 153:    if statement
  set_difference(allNodes.begin(),allNodes.end(),
		 allAssignedProbNodes.begin(),allAssignedProbNodes.end(),
		 inserter(nodesThatGiveNoProb,
			  nodesThatGiveNoProb.end()));

  if (nodesThatGiveNoProb.size() > 0) {
    fprintf(stderr,"INTERNAL ERROR: some nodes do not give probability to any clique, those nodes include: ");
    printRVSet(stderr,nodesThatGiveNoProb);
    coredump("Possibly corrupt trifile. Exiting program ...");
  }
#endif
}

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::assignRVsToCliques()
 *    Assign to the cliques the random variables that will
 *    be used to generate the tables for each clique, and do
 *    it for a particular section 'section' starting at
 *    the clique number given by rootClique.
 * 
 *    The assumption is that collectEvidence will be called with this
 *    rootClique, so the assignment of variables to cliques is an attempt
 *    to assign parents before children, and assign parents farther away
 *    from the root than children, so that by the time children are encountered,
 *    the parent values are already enumerated.
 *
 *    The basic algorithm is recursive, and uses a helper routine for
 *    that purpose. It descends down the graph in depth first order,
 *    attempting to assign a random variable to a clique. Note that
 *    it is possible for a rv to not be assigned (i.e., if no clique
 *    exist that contains all the rvs parents). If this happens, it
 *    is not an error, rather it should be the case that the rv lives with
 *    its parents in another section.
 *
 * Preconditions:
 *     createSectionJunctionTrees(), computeSectionInterfaces() must
 *     have been called. The section must not be empty.
 *
 * Postconditions:
 *     the section so specified in the code below
 *     are now such that their cliques have their 'assignedNodes' 
 *     member filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::assignRVsToCliques(const char *const sectionName,
				     JT_Partition &section,
				     const unsigned rootClique,
				     const char* varSectionAssignmentPrior,
				     const char *varCliqueAssignmentPrior)
{
  vector<RV*> sortedNodes;

  // We use a constrained topological sort based on command line
  // arguments, so that variables are considered in what hopefully
  // will be a good order (e.g., might have it such that continuous
  // variables or discrete observations come as early as possible in
  // the ordering).  This will also be beneficial when sampling hidden
  // continuous nodes.
  infoMsg(IM::Giga,"Sorting section variables using priority order (%s)\n",
	  (varSectionAssignmentPrior?varSectionAssignmentPrior:"NULL"));

  GraphicalModel::topologicalSortWPriority(section.nodes,section.nodes,sortedNodes,varSectionAssignmentPrior);

  // printf("have %d sorted nodes and %d cliques\n",sortedNodes.size(),section.cliques.size());

  // update the cumulative RV assignments before we begin.
  getCumulativeAssignedNodes(section,rootClique);

  for (unsigned n=0;n<sortedNodes.size();n++) {

    RV* rv = sortedNodes[n];

    // precompute the parents set for set intersection operations.
    // The one stored in the rv is a vector, but we need a set, so we
    // pre-compute it here. While we're at it, we compute if all
    // parents are observed since the behaviour is a bit
    // different in this case.
    // TODO: check if this is really necessary, as in template routines vector is a container
    //       so can be used as needed.
    set<RV*> parSet;
    bool allParentsObserved=true;
    for (unsigned p=0;p<rv->allParents.size();p++) {
      // TODO: see if there is a faster STL way to go from vector to set.
      // fprintf(stderr,"about to insert rv with address %X\n",(void*)rv->allPossibleParents[p]);
      parSet.insert(rv->allParents[p]);
      if (rv->allParents[p]->hidden())
	allParentsObserved = false;
    }

    // create a map to keep track of the scores of a given clique in
    // which RV should produce probability.  This utilizes the fact
    // that the multimap stores values in *ascending* order based on
    // key (here the score), and so the first one (i.e.,
    // scoreSet.begin()) should have the lowest weight.
    multimap < vector<double>, unsigned> scoreSet;

    unsigned numberOfTimesAssigned = 0;
    assignRVToClique(sectionName,
		     section,rootClique,
		     0,
		     rv,
		     numberOfTimesAssigned,
		     parSet,allParentsObserved,
		     scoreSet);
    infoMsg(IM::Max,
	    "Section %s: random variable %s(%d) with its parents contained in %d cliques\n",
	    sectionName,
	    rv->name().c_str(),rv->frame(),numberOfTimesAssigned);

    if (numberOfTimesAssigned > 0) {
      // then it was assigned in this section at least once, but we
      // still need to check if it should give probability to some
      // clique within this section. We need to check both if the
      // scoreSet has a value, and also if it was not assigned, since
      // the scoreSet might have obtained a value from a branch of the
      // JT below which didn't have this rv assigned with probability.

      if ((scoreSet.size() > 0)  &&
	  (section.cliques[rootClique].cumulativeAssignedProbNodes.find(rv) == 
	   section.cliques[rootClique].cumulativeAssignedProbNodes.end())) {

	// Then the node will be a probability contributer to one clique
	// in this section.

	// We choose one of those cliques to be the one the node
	// contributes probabilty to.
	unsigned clique_num = (*(scoreSet.begin())).second;
	section.cliques[clique_num].assignedProbNodes.insert(rv);
	if (IM::messageGlb(IM::Max)) {
	  infoMsg(IM::Max,
		  "Section %s: random variable %s(%d) giving probability to clique %d.\n",
		  sectionName,
		  rv->name().c_str(),rv->frame(),clique_num);
	}
      }

      // update the cumulative RV assignments.
      getCumulativeAssignedNodes(section,rootClique);

    } else {

      // rv was not assigned to this section, it must be the case
      // that it will be assigned to a different section. This
      // could come from the section before or the section after
      // the current section. If there are only forward-directed
      // arrows, it will come from the previous section (and vice
      // versa if there are only backwards going edges). If we have
      // both directed edges, it could be assigned in either left or
      // right section.

      // Note that it is impossible for a node to be assigned to give
      // probabiltiy in two sections. The reason is the
      // following. The only nodes that live in both sections are
      // the interface nodes.  These nodes have been special cased
      // above using the 'alreadyAProbContributer' variable, and
      // since we keep cumulative track of nodes that are to
      // provide probability, we won't do this more than once
      // across sections.

      infoMsg(IM::Max,"Section %s: random variable %s(%d) not assigned in current section\n",sectionName,
	      rv->name().c_str(),rv->frame());
      // in any event, keep track of this node, perhaps useful for jtWeight.
      section.unassignedInPartition.insert(rv);
    }
  }
}




/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::assignRVToClique()
 *    Helper routine for assignRVsToCliques() above.
 *
 * Preconditions:
 *     Must be called only by assignRVsToCliques().
 *
 * Postconditions:
 *     the random variable might be assigned (but the assignment
 *     might also be done by the caller).
 *
 * Side Effects:
 *     Might change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::assignRVToClique(const char *const sectionName,
				   JT_Partition& section,
				   const unsigned root,
				   const unsigned depth,
				   RV* rv,
				   unsigned& numberOfTimesAssigned,
				   set<RV*>& parSet,
				   const bool allParentsObserved,
				   multimap < vector<double>, unsigned >& scoreSet)
{
  // keep a reference for easy access
  MaxClique& curClique = section.cliques[root];

  // First, check if directed_closure(rv) == (rv U parents(rv)) not in
  // the current clique, and if so then continue on down.
  bool closure_in_clique = (curClique.nodes.find(rv) != curClique.nodes.end());

  // printf("clique%d: from rv, closure_in_clique = %d\n",root,closure_in_clique);

  if (closure_in_clique && !allParentsObserved) {
    // need to further check that parents are in clique
    set<RV*> parentsInClique;
    set_intersection(curClique.nodes.begin(),
		     curClique.nodes.end(),
		     parSet.begin(),parSet.end(),
		     inserter(parentsInClique,parentsInClique.end()));

    // First check for equal size. If the size is equal, then *all*
    // parents are in the clique.
    closure_in_clique = (parentsInClique.size() == parSet.size());
    if (!closure_in_clique) {
      // rather than only checking for equal size, we here also check
      // if no *hidden* parents are out of the clique, since that is
      // all we need to assign this node to this clique. In other
      // words, if all parents that are out of clique are observed, we
      // can still go. In fact, assuming moralization doesn't moralize
      // w.r.t. an observed parent, this step is necessary sometimes
      // to get the child assigned.

      assert ( parentsInClique.size() < parSet.size() );

      set<RV*> parentsOutOfClique;
      set_difference(parSet.begin(),parSet.end(),
		     parentsInClique.begin(),parentsInClique.end(),
		     inserter(parentsOutOfClique,parentsOutOfClique.end()));

      // Check to see if all parents out of the clique are observed,
      // and if so, the closure is in the clique.
      set<RV*>::iterator it;
      closure_in_clique = true;
      for (it = parentsOutOfClique.begin(); it != parentsOutOfClique.end(); it++) {
	RV* rv = (*it);
	if (rv->hidden()) {
	  closure_in_clique = false;
	  break;
	}
      }
    }
  }

  // printf("clique%d: from par(rv), closure_not_in_clique = %d, num
  // children =
  // %d\n",root,closure_not_in_clique,section.cliques[root].children.size());

  if (closure_in_clique) {
    // So closure(rv) is in current clique, meaning rv and *all* its
    // hidden parents are members of this clique, but it is not
    // nec. the case that all of rv's parents are themselves assigned
    // to produce probabilities in this clique.

    // Since closure is in the clique, we "assign" this node to this
    // clique. Note, this does not necessarily mean that the node will
    // contribute probabilties to this clique's potential function,
    // only that we will use the RV iterator in some cases (sparse,
    // deterministic, sampling) to iterate over this RV after its
    // parents have been assigned (as long as the RV doesn't come in
    // as a separator that is).
    curClique.assignedNodes.insert(rv);
    // we include this in the sorted assigned nodes for now, even though
    // it probably will be re-sorted at a later time.
    curClique.sortedAssignedNodes.push_back(rv);
    numberOfTimesAssigned++;

    infoMsg(IM::Max,
	    "Section %s: RV %s(%d) and its parents in clique %d\n",
	    sectionName,
	    rv->name().c_str(),rv->frame(),root);

    // We don't make the probability assignment decision until after
    // we consider all possible such candidate cliques.  But first,
    // check if already set up to give probability to some other
    // clique.
    if (curClique.cumulativeAssignedProbNodes.find(rv) == curClique.cumulativeAssignedProbNodes.end()) {
      // So RV is not yet giving probability anywhere. We deal with
      // the score stuff to find a good clique for this variable to
      // assign probability.

      infoMsg(IM::Max,
	      "Section %s: considering clique %d for random variable %s(%d) probability contribution\n",
	      sectionName,
	      root,
	      rv->name().c_str(),rv->frame());

      // In this discussion that follows, we make a distinction
      // between node parents and children (i.e., original graph
      // parents and children) and JT parents and children which are
      // cliques in the JT. A child clique c of a parent p in the JT
      // is one such that c ~ p and that depth(c) = depth(p)+1, where
      // depth(root)=0 in the junction tree.

      // The goal here is to figure out in which of a set of possible
      // cliques the node should be producing probabilities.

      // TODO: make these various options available in priority order
      // to the command line, just like the topologically-constrainted
      // variable within-clique iteration order, currently available
      // via -vcap.


      // Since the multimap stores values in *ascending* order based
      // on key, lower numbers are prefered.

      // Compute (and ultimately push back) items in decreasing order of
      // priority.  Lower numbers is better (e.g., more negative or less
      // positive is higher priority).  First thing inserted has highest
      // priority.  At some time later, the rv will be assigned to the
      // clique with the *LOWEST* score.

      // What we do, is continue on down. It may be the the case that we
      // find a clique with rv and all of rv's node parents assigned (in
      // which case we assign rv right then), but in the mean time we
      // keep track of this clique's "score" using a set of
      // heuristics. 

      // TODO: change it so that only heuristics that are
      // used are computed.

      // Compute number of previous parents with their probability assigned.
      // The heuristic here is, have this node contribute probabiltiy
      // to this clique if many of its parents are already doing so, which
      // might produce a clique with good pruning behavior.
      count_iterator<set <RV*> > res;
      set_intersection(curClique.assignedProbNodes.begin(),
		       curClique.assignedProbNodes.end(),
		       parSet.begin(),parSet.end(),
		       res);
      int num_parents_with_probability = (int) res.count();
      // Previous Parents (earlier in the junction tree) with their
      // probabilities in Junction Tree.  We add this to the above.
      set_intersection(curClique.cumulativeAssignedProbNodes.begin(),
		       curClique.cumulativeAssignedProbNodes.end(),
		       parSet.begin(),parSet.end(),
		       res);
      num_parents_with_probability += (int) res.count();
      // negate so that lower is preferable.
      num_parents_with_probability *= -1;

      // Distance from root, among the higher priorities that are equal,
      // try to be as far away from the root as possible so as to prune
      // away as much zero as possible as early as possible. I.e., try
      // encourage it to have as much "influence" as possible in JT
      // parents.
      // negate so that lower (farther from the root) is preferable.
      const int distance_from_root = -(int)depth;

      // Number of children in current clique. If rv has lots of
      // children in this clique, it is hopeful that other parents of
      // those children might also be assigned to the same clique.
      int numChildren = 0;
      for (unsigned i=0;i<rv->allChildren.size();i++) {
	RV* child = rv->allChildren[i];
	if (curClique.nodes.find(child) != curClique.nodes.end())
	  numChildren++;
      }
      // negate so that lower (more children in current clique) is preferable.
      numChildren *= -1;

      // when the node is continuous, assign it to a clique that is the
      // smallest possible in terms of weight.
      // Again, negate so that larger weight is preferable.
      double weight = - curClique.weight();

      // TODO: add heuristic that assigns RV to smallest clique, in
      // terms of size. This might be useful for EM training,
      // particularly of continuous Gaussian variables, in that only
      // the parents of the Gaussian will be iterated over, rather
      // than parents of parents.

      // And so on. We can create as many heuristics as we want.
      // alternatively, perhaps take a weighted average of some of them.
      // TODO: add more heuristics here, and/or produce better
      // prioritized order above.

      // Other possible heuristics:
      //    1) a clique where it has the greatest number of node
      //    descendants (but this is only a heuristic since to gain
      //    benefit, we would need to have it be such that those node
      //    descendants have their node parents in clique as well.


      // We've now computed a bunch of heursitcs, push them
      // into an array in priority order to be scored later
      vector<double> score;

      // most important, distance from root, fail first, prune first.
      score.push_back(distance_from_root);
      // next most, break ties with weight of clique, encourage fail first principle
      score.push_back(weight);
      // next, if lots of parents in clique, adding this here might help prune.
      score.push_back(num_parents_with_probability);      
      // break ties with number of children, likely to help others prune.
      score.push_back(numChildren);

#if 0
      if (rv->discrete()) {
	DiscRV* drv = (DiscRV*)rv;
	// 1st thing pushed has highest priority.
	if (drv->sparse() || drv->deterministic()) {
	  // if the RV is sparse/deterministic, use distance from root
	  // as the higest priority, since a sparse node will have
	  // lots of zeros, regardless of beam width. Goal is to prune
	  // away as early as possible.
	  score.push_back(distance_from_root);
	  score.push_back(num_parents_with_probability);
	} else {
	  // if the RV is not sparse, try to assign it to a clique
	  // with as many of its parents assigning probabiltiy, since
	  // a it will make for a clique that will prune more
	  // effectively
	  score.push_back(num_parents_with_probability);
	  score.push_back(distance_from_root);
	}
	// Last, assign based on prefering cliques where this variable
	// has lots of children (the children might then use the fact
	// that this RV gives probability to encourage it to live here
	// as well, again to improve pruning performance).
	score.push_back(numChildren);
      } else {
	// RV is continue. Pefer small weight.
	score.push_back(weight);
	// RV is continue. Then perfer lots of parents with probability.
	score.push_back(num_parents_with_probability);      
      }
#endif

      // done inserting heuristicss, now insert the score and the
      // current clique into the set.
      pair < vector<double>, unsigned> p(score,root);
      scoreSet.insert(p);
    }
  }

  // continue on down.
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    const unsigned child = section.cliques[root].children[childNo];
    assignRVToClique(sectionName,section,child,depth+1,rv,numberOfTimesAssigned,parSet,allParentsObserved,scoreSet);
  }

}



/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::computeUnassignedCliqueNodes()
 *    Go through each clique and compute the unassigned nodes.
 *    This is only done when ceSeparatorDrivenInference == false.
 *
 * Preconditions:
 *     Must be called in appropriate order (setUpDataStructures).
 *
 * Postconditions:
 *     Cliques have this member variable filled in.
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::computeUnassignedCliqueNodes(JT_Partition& part)
{
  if (true /* MaxClique::ceSeparatorDrivenInference == false */ ) {
    // these are only needed by clique driven inference.
    for (unsigned cliqueNum=0;cliqueNum<part.cliques.size();cliqueNum++) {
      part.cliques[cliqueNum].computeUnassignedCliqueNodes();
    }
  }
}
void
SectionScheduler::computeUnassignedCliqueNodes()
{
  computeUnassignedCliqueNodes(P1);
  computeUnassignedCliqueNodes(Co);
  computeUnassignedCliqueNodes(E1);
}





/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::assignFactorsCliques()
 *    Assign any existing factors to cliques in the junction tree.
 *
 * Preconditions:
 *     createSectionJunctionTrees(), computeSectionInterfaces(), and
 *     createDirectedGraphOfCliques() must have been called before 
 *     this can be called.
 *
 * Postconditions:
 *     each of the factors (if any) in the current section
 *     have now been assigned to the clique that is guaranteed
 *     to contain that factor. There must exist a clique that
 *     contains the factor by the triangulation code (there'll be
 *     an error if we can't find a clique to place the factor).
 *
 * Side Effects:
 *     Will change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::assignFactorsToCliques() {
  assignFactorsToCliques(P1);
  assignFactorsToCliques(Co);
  assignFactorsToCliques(E1);
}

void 
SectionScheduler::assignFactorsToCliques(JT_Partition& section) {
  // for each factor, assign it to any clique that contains
  // it.

  for (unsigned f=0;f<section.factorCliques.size();f++) {
    // temporary: for now, only assign factors that
    // correspond ot the allequal constriant.
    if ((section.factorCliques[f].factorInfo->fType != FactorInfo::ft_symmetricConstraint)
	||
	(section.factorCliques[f].factorInfo->symmetricConstraintInfo.symmetricConstraintType
	 != FactorInfo::sct_allVarsEqual))
      continue;

    for (unsigned c=0;c<section.cliques.size();c++) {
      infoMsg(IM::Giga,"Checking if factor %d fits in clique %d, in section\n",
	      f,c);
      if (firstRVSetContainedInSecond(section.factorCliques[f].nodes,
				      section.cliques[c].nodes)) 
      {
	// then this factor becomes assigned to this cliuqe.
	section.cliques[c].assignedFactors.push_back(f);
	infoMsg(IM::Giga,"Assigning factor %s(%d) to clique %d\n",
		section.factorCliques[f].factorInfo->name.c_str(),
		section.factorCliques[f].factorInfo->frame,
		c);
      }
    }
  }
}




/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::getCumulativeAssignedNodes()
 *    Decends down the directed JT graph and computes and
 *    assignes the variable cumulativeAsignedNodes.
 *
 * Preconditions:
 *     Must be called only by assignRVsToCliques()
 *
 * Postconditions:
 *     the clique's member variable cumulativeAsignedNodes has
 *     been re-computed.
 *
 * Side Effects:
 *     Might change the clique data member variables
 *
 * Results:
 *     None.
 *
 *-----------------------------------------------------------------------
 */

void
SectionScheduler::getCumulativeAssignedNodes(JT_Partition &section, vector<unsigned> const &roots) {
  for (unsigned i=0; i < roots.size(); ++i) {
    getCumulativeAssignedNodes(section, roots[i]);
  }
}
void
SectionScheduler::getCumulativeAssignedNodes(JT_Partition& section,
					     const unsigned root)
{
  MaxClique& curClique = section.cliques[root];

  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = curClique.children[childNo];
    getCumulativeAssignedNodes(section,child);
    // Note: this will do a bunch of redundant work (i.e., inserting
    // elements into sets where the elements are already contained in
    // the sets), but it doesn't need to run that fast, so the
    // redundant work is not a big deal.
    unionRVs(section.cliques[child].cumulativeAssignedNodes,
	     section.cliques[child].assignedNodes,
	     curClique.cumulativeAssignedNodes,true);
    unionRVs(section.cliques[child].cumulativeAssignedProbNodes,
	     section.cliques[child].assignedProbNodes,
	     curClique.cumulativeAssignedProbNodes,true);

  }
}



/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::getCumulativeUnassignedIteratedNodes()
 *
 *   Computes for each clique the union of the set of nodes which have
 *   been iterated unassigned in previous cliques in the JT. This is
 *   used to determine which, in each clique, nodes should be iterated
 *   over and which are already assigned by a separator driven
 *   iteration.
 *
 * Preconditions:
 *   computeSeparatorIterationOrder() must have been called. Section
 *   in argument 'section' must not be empty.
 *
 * Postconditions:
 *   cliques know which nodes to iterate over or not.
 *
 * Side Effects:
 *    Changes member variables within cliques of this section.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::getCumulativeUnassignedIteratedNodes(JT_Partition &section,
						       const unsigned root)
{
  MaxClique& curClique = section.cliques[root];
  set<RV*>& res = curClique.cumulativeUnassignedIteratedNodes;
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {

    const unsigned child = curClique.children[childNo];
    getCumulativeUnassignedIteratedNodes(section,child);
    set_union(section.cliques[child].cumulativeUnassignedIteratedNodes.begin(),
	      section.cliques[child].cumulativeUnassignedIteratedNodes.end(),
	      section.cliques[child].unassignedIteratedNodes.begin(),
	      section.cliques[child].unassignedIteratedNodes.end(),
	      inserter(res,res.end()));
  }
  // curClique.computeAssignedNodesDispositions();
  // if there have been any changes, then we need to re-sort
}


// Computes the preceding iterated unassigned nodes and therein the
// set of assigned nodes in each clique that should/shouldn't be
// iterated.
void
SectionScheduler::getCumulativeUnassignedIteratedNodes()
{
#if 0
  set<RV*> res;

  if (P1.cliques.size() > 0) {
    // TODO: this should be a member function in P1.
    getCumulativeUnassignedIteratedNodes(P1,P_ri_to_C);
    res.clear();
    set_union(P1.cliques[P_ri_to_C].unassignedIteratedNodes.begin(),
	      P1.cliques[P_ri_to_C].unassignedIteratedNodes.end(),
	      P1.cliques[P_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),
	      P1.cliques[P_ri_to_C].cumulativeUnassignedIteratedNodes.end(),	    
	      inserter(res,res.end()));
  }

  // Co is never empty.
  Co.cliques[C_li_to_P].cumulativeUnassignedIteratedNodes = res;
  getCumulativeUnassignedIteratedNodes(Co,C_ri_to_C);

  if (E1.cliques.size() > 0) {
    res.clear();
    set_union(Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].unassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),
	      Co.cliques[C_ri_to_C].cumulativeUnassignedIteratedNodes.begin(),	    
	      inserter(res,res.end()));
    E1.cliques[E_li_to_C].cumulativeUnassignedIteratedNodes = res;
    getCumulativeUnassignedIteratedNodes(E1,E_root_clique);
  }
#endif
}




// For the three sections, set up the different message passing
// orders that are to be used. This basically just does a tree
// traversal using the previously selected root.

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::setUpMessagePassingOrders()
 *   sets up the message passing orders for the various
 *   basic sections.
 * 
 *
 * Preconditions:
 *   Sections must be instantiated, and interface cliques
 *   must have been computed (i.e., computeSectionInterfaces()
 *   must have been called)
 *
 * Postconditions:
 *   The member variables P_to_C_message_order (etc.) are created
 *   with the order of message passing.
 *
 * Side Effects:
 *    Changes member variables as described above.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::setUpMessagePassingOrders() {
#if 0
  setUpMessagePassingOrder(P1,
			   P_ri_to_C,
			   P1_message_order,
			   ~0x0,
			   P1_leaf_cliques);
  setUpMessagePassingOrder(Co,
			   C_ri_to_C,
			   Co_message_order,
			   C_li_to_P,
			   Co_leaf_cliques);
  setUpMessagePassingOrder(E1,
			   E_root_clique,
			   E1_message_order,
			   E_li_to_C,
			   E1_leaf_cliques);
#endif
}

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::setUpMessagePassingOrders()
 *   sets up the message passing orders for a given
 *   specified section, and places the result in the argument.
 *   Roots the undirected clique tree at 'root' thereby creating
 *   a directed tree of cliques. This does a recursive DFS tree traversal
 *   start at root to make the order, and places the visitation
 *   order in the array.
 * 
 * Preconditions:
 *   Same as setUpMessagePassingOrders() but for the
 *   particular arguments being passed in.
 *
 * Postconditions:
 *   order is set up.
 *.
 * Side Effects:
 *    none
 *
 * Results:
 *    None, but passed back in arguments.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::setUpMessagePassingOrder(JT_Partition &section,
					   const unsigned root,
					   vector< pair<unsigned,unsigned> >&order,
					   const unsigned excludeFromLeafCliques,
					   vector< unsigned>& leaf_cliques) 
{
  // first, handle case for potential empty P or E section. 
  if (section.cliques.size() == 0)
    return;

  order.clear();
  // a tree of N nodes has N-1 edges.
  order.reserve(section.cliques.size() - 1);
  setUpMessagePassingOrderRecurse(section,root,order,excludeFromLeafCliques,leaf_cliques);
}

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::setUpMessagePassingOrderRecurse()
 *
 *   Recursive support routine for setUpMessagePassingOrder() and
 *   should only be called from there. See that
 *   routine for documentation.
 * 
 *-----------------------------------------------------------------------
 */					
void
SectionScheduler::setUpMessagePassingOrderRecurse(JT_Partition &section,
					      const unsigned root,
					      vector< pair<unsigned,unsigned> >& order,
					      const unsigned excludeFromLeafCliques,
					      vector < unsigned >& leaf_cliques)
{
  if (section.cliques[root].children.size() == 0) {
    // store leaf nodes here if children.size() == 0.
    if (root != excludeFromLeafCliques)
      leaf_cliques.push_back(root);
  } else {
    for (unsigned childNo=0;
	 childNo<section.cliques[root].children.size();
	 childNo++) {
      const unsigned child = section.cliques[root].children[childNo];
      setUpMessagePassingOrderRecurse(section,child,order,
				      excludeFromLeafCliques,leaf_cliques);
      pair<unsigned,unsigned> msg;
      msg.first = child;
      msg.second = root;
      order.push_back(msg);
    }
  }
}


// Separator creation, meaning create the seperator objects
// both within and between sections. Given two neighboring
// sections L and R, the separator between the interface
// cliques in L and R is contained in R.


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createSeparator()
 *   Creates the separator object based on the message passing order
 *   that has been created at this point for this section.
 *
 * Preconditions:
 *   setUpMessagePassingOrder() must have already been called for
 *   this section.
 *
 * Postconditions:
 *   The separators for this section have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the sections separators data structures.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::createSeparators(JT_Partition &section, vector< pair<unsigned,unsigned> > &order) {
  section.separators.clear();

  // first set up all the cliques so that they know who their separators are.
  for (unsigned p=0; p < order.size(); p++) {
    const unsigned from = order[p].first;
    const unsigned to   = order[p].second;    
    const unsigned sepNo = section.separators.size();
    section.separators.push_back(SeparatorClique(section.cliques[from],section.cliques[to]));
    section.cliques[from].ceSendSeparators.push_back(sepNo);
    section.cliques[to].ceReceiveSeparators.push_back(sepNo);
  }
}

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createSeperators()
 *
 *   Creates the seperator objects based on the message passing order
 *   that has been created at this point. It creates seperators for
 *   the base sections P1, Co, and E1. Each seperator contains nodes
 *   consisting of the intersection of the maxCliques between which a
 *   message is sent.
 *
 *   The order the separators in the separator array are created in each section are:
 *      1) first come the separators between cliques entirely contained within
 *         the section.
 *      2) next come any VE separators, if there are any.
 *      3) last comes the left-interface (LI) separator, which is the
 *         separator between the right-interface clique in the left
 *         section and the left-interface clique in the current
 *         section. Note that P1 does not have such an entry
 *         in its separator array.
 *   
 *   See comment in class JT_Partition in .h file at the 'separators'
 *   member for more information on separators and the order that
 *   they come in.
 * 
 * Preconditions:
 *   setUpMessagePassingOrder() must have already been called.
 *
 * Postconditions:
 *   The seperators for sections[0-3] have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the sections seperators data structures.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::createSeparators() {
  if (P1.cliques.size() > 0) { 
    createSeparators(P1, P1_message_order);

    // set to invalid value, since the ri clique sends to a separator
    // in the next section, rather than internal to this section.

    // placeholders for outgoing separators from P to C
    for (unsigned i=0; i < P1.section_ri.size(); ++i) {
      P1.cliques[ P1.section_ri[i] /* P_ri_to_C */ ].ceSendSeparators.push_back(~0x0); 
    }

    // insert VE seps last
    if (useVESeparators && (veSeparatorWhere & VESEP_WHERE_P)) {
      infoMsg(IM::Max,"Searching for VE seps in P1\n");
      createVESeparators(P1);
    } else
      P1.numVEseps = 0;

    // P1 has no LI, so nothing more to do for P1 here.

  }

  // Cs are never empty.
  assert ( Co.cliques.size() > 0 );

  // first create the normal seprators
  createSeparators(Co, Co_message_order);

  // Next, potentially insert VE seps, before the last one.
  if (useVESeparators && (veSeparatorWhere & VESEP_WHERE_C)) {
    infoMsg(IM::Max,"Searching for VE seps in Co\n");
    createVESeparators(Co);
  } else
    Co.numVEseps = 0;

  // Final separator is guaranteed *always* to be the interface
  // separator to the left section, the LI separator.
  if (gm_template.usesLeftInterface()) {
    // left interface case
    assert ( E1.cliques.size() > 0 );

    if (P1.cliques.size() > 0) {
      // Create separator of interface cliques, specifically the
      // separator between, P1's right-interface clique and Co's left
      // interface clqiue. 
      // 
      // Note that the interface separator always has to be the last
      // separator inserted.

#if 0
// FIXME - This will be problematic for supporting multiple inference
//         architectures at inference time: ia_name:C knows which clique #s
//         to use for its interfaces, but it won't know anything about the
//         cliques in the sections to the left or right being interfaced to/from.
//         Could either dynamically construct the separators at inference time
//         once the other section's inference architecture is known (probably
//         expensive), or pre-compute all possible X:P_to_Y:C and Y:C_to_Z:E
//         combinations at use the proper one once the inference architecture
//         being interfaced with is known.

      Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_to_C],Co.cliques[C_li_to_P]));
      // update right sections LI clique to include new separator
      Co.cliques[C_li_to_P].ceReceiveSeparators.push_back(Co.separators.size()-1);
#else

      if (P1.section_ri.size() != Co.section_li.size()) {
	// FIXME - add ia_names to error message
	error("ERROR: PC interface size mismatch\n");
      }
      
      // build PC separators in C (i.e., P's right & C's left interface separators)
      for (unsigned i=0; i < P1.section_ri.size(); ++i) {
	unsigned P_ri_clique_num = P1.section_ri[i], 
                 C_li_clique_num = Co.section_li[i];
	Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_clique_num], 
						Co.cliques[C_li_clique_num]));
	Co.cliques[C_li_clique_num].ceReceiveSeparators.push_back(Co.separators.size()-1);
      }
#endif
    } else {
      // Then P1 is empty. We insert a dummy separator that has
      // been set up so that C'_i can go to C'_{i+1}. For i>0 this
      // works fine, for i=0 we will need to make sure during
      // inference not to use this separator since it will have
      // nothing in it.

      if (Co.section_li.size() > 0) {
	// then there will at least be a Ci <-> C_{i+1}
	// intersection. I.e., we are guaranteed that the chunks are
	// not disconnected from each other. This can happen if C
	// reaches into the future (i.e., reaches into its next self
	// and into E) --- for example, if there are children
	// variables in Ci with their parents in C_{i+1}. It is not
	// possible in this case for a child in Ci to have a parent in
	// C_{i-1} since we're in the size(P) == 0 case here.
	for (unsigned i=0; i < Co.section_li.size(); ++i) {
	  // FIXME - the GMTK 1.X code has C_li_to_P -> C_li_to_P ?
	  unsigned C_ri_clique_num = Co.section_ri[i],
                   C_li_clique_num = Co.section_li[i];
	  Co.separators.push_back(SeparatorClique(Co.cliques[C_ri_clique_num],
						  Co.cliques[C_li_clique_num]));
	  Co.cliques[C_li_clique_num].ceReceiveSeparators.push_back(Co.separators.size()-1);
	}
      } else {
	// then there will not be a Ci <-> C_{i+1} intersection.  We
	// need to create an empty separator since the chunks are
	// disconnected from each other.
	set<RV*> empty;
	MaxClique mc(empty);
	Co.separators.push_back(SeparatorClique(mc,mc));
	for (unsigned i=0; i < Co.section_li.size(); ++i) {
	  // update right sections LI clique to include new separator
	  Co.cliques[Co.section_li[i]].ceReceiveSeparators.push_back(Co.separators.size()-1);
	}
      }
    }
  } else {
    // right interface case
    assert ( P1.cliques.size() > 0 );

    if (P1.section_ri.size() != Co.section_li.size()) {
      // FIXME - add ia_names to error message
      error("ERROR: PC interface size mismatch\n");
    }

    // Note that the interface separator always has to be the last separator inserted.
    for (unsigned i=0; i < Co.section_li.size(); ++i) {
      unsigned P_ri_clique_num = P1.section_ri[i],
	       C_li_clique_num = Co.section_li[i];
      Co.separators.push_back(SeparatorClique(P1.cliques[P_ri_clique_num],
					      Co.cliques[C_li_clique_num]));
      // update right sections LI clique to include new separator
      Co.cliques[C_li_clique_num].ceReceiveSeparators.push_back(Co.separators.size()-1);
    }
  }

  // don't update left sections RI clique's send separator since handled explicitly
  for (unsigned i=0; i < Co.section_ri.size(); ++i) {
    Co.cliques[ Co.section_ri[i] ].ceSendSeparators.push_back(~0x0); //set to invalid value
  }

  if (E1.cliques.size() > 0) { 

    // normal separators
    createSeparators(E1, E1_message_order);


    // Next, potentially insert any VE seps, before final one.
    if (useVESeparators  && (veSeparatorWhere & VESEP_WHERE_E)) {
      infoMsg(IM::Max,"Searching for VE seps in E1\n");
      createVESeparators(E1);
    } else
      E1.numVEseps = 0;

    // Final LI separator one is guaranteed *always* to be the
    // interface separator to the left section.

    if(Co.section_ri.size() != E1.section_li.size()) {
      // FIXME - add ia_names to error message
      error("ERROR: CE interface size mismatch\n");
    }
    // 
    // Create separator of interface cliques. Co is never empty. Note that the interface
    // separator always has to be the last separator inserted.
    for (unsigned i=0; i < Co.section_ri.size(); ++i) {
      unsigned C_ri_clique_num = Co.section_ri[i],
	       E_li_clique_num = E1.section_li[i];
      E1.separators.push_back(SeparatorClique(Co.cliques[C_ri_clique_num],
					      E1.cliques[E_li_clique_num]));
      // update right sections LI clique to include new separator
      E1.cliques[E_li_clique_num].ceReceiveSeparators.push_back(E1.separators.size()-1);
    }
    // don't update left sections RI clique's send separator since handeled explicitly
    for (unsigned i=0; i < E_root_clique.size(); ++i) {
      E1.cliques[E_root_clique[i]].ceSendSeparators.push_back(~0x0); //set to invalid value
    }
  }
}


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::createVESeparator()
 *   Creates the virtual evidence (VE) separator objects based on any
 *   virtual evidence that might live in the current cliques.
 *
 * Preconditions:
 *   createSeparators() must have already been called.
 *
 * Postconditions:
 *   The VE separators for this section have been created
 *   the right and left.
 *
 * Side Effects:
 *   Modifies the sections separators data structures. Pushes
 *   each such separator back at the end of separator list.
 *
 * Results:
 *   none
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::createVESeparators(JT_Partition& section) {
  section.numVEseps = 0;
  if (!useVESeparators) {
    return;
  }
  for (unsigned c=0;c<section.cliques.size();c++) {
    // there might several VE separators.
    unsigned numVEsepsInClique = section.cliques[c].computeVESeparators();
    if (numVEsepsInClique > 0) {
      if (IM::messageGlb(IM::Max)) {
	printf("Found %d VE seps in clique number %d of current section\n",numVEsepsInClique,c);
	for (unsigned i=0;i<numVEsepsInClique;i++) {
	  printf("VE sep %d,",i);
	  if (section.cliques[c].veSeparators[i].grandChild != NULL) 
	    printf("pcg:");
	  else
	    printf("pc:");
	  printRVSetAndCards(stdout,section.cliques[c].veSeparators[i].parents,false);
	  printf(",c=");
	  RV2DRV(section.cliques[c].veSeparators[i].child)->printNameFrameCard(stdout,false);
	  if (section.cliques[c].veSeparators[i].grandChild != NULL) {
	    printf(",gc=");
	    RV2DRV(section.cliques[c].veSeparators[i].grandChild)->printNameFrameCard(stdout,false);
	  }
	  printf("\n");
	}
      }
    }
    for (unsigned n=0;n<numVEsepsInClique;n++) {
      const unsigned sepNo = section.separators.size();
      section.separators.push_back(SeparatorClique(section.cliques[c].VESeparatorInfo(n)));
      section.cliques[c].ceReceiveSeparators.push_back(sepNo);
    }
    section.numVEseps += numVEsepsInClique;
  }
}


// Separator iteration order and accumulated set intersection
// creation for separator driven clique potential creation, and
// also updates the seperators partial accumulator structure and
// sets up cliques other variables.

/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::computeSeparatorIterationOrders()/computeSeparatorIterationOrder()
 *
 *   1) computes the order in which we iterate over the separators
 *   when we do seperator driven clique potential function creation,
 *   for the given section.
 *
 *   2) Also, set up the separator structures so they know the
 *   intersection of each of their nodes with the accumulated union of
 *   all nodes in separators earlier in the order.
 *
 *   3) Also, assigns unassignedIteratedNodes in this clique.
 *
 * Preconditions:
 *   separators must be created, should only be called from
 *   the general routine above w/o arguments.
 *
 * Postconditions:
 *   Separator orders have now been created within given section.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques of each section.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::computeSeparatorIterationOrders(JT_Partition& section) {
  for (unsigned cliqueNum=0; cliqueNum < section.cliques.size(); cliqueNum++) {
    computeSeparatorIterationOrder(section.cliques[cliqueNum],section);
  }
}

void 
SectionScheduler::computeSeparatorIterationOrder(MaxClique& clique, JT_Partition& section) {
  // do an order based on the MST over the separators.
  const unsigned numSeparators = clique.ceReceiveSeparators.size();

  // Build up the partial and then ultimately the final union of of
  // all nodes in separators for incomming messages for this clique.
  clique.unionIncommingCESeps.clear();

  if (numSeparators == 0) {
    // This must be a leaf-node clique relatve to root.
    // 'unionIncommingCESeps' is already empty so no need to do anything there.
  } else if (numSeparators == 1) {
    // shortcut to separator 0
    SeparatorClique& s0 = section.separators[clique.ceReceiveSeparators[0]];
    clique.unionIncommingCESeps = s0.nodes;
    s0.accumulatedIntersection.clear();
    s0.remainder = s0.nodes;
    assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
  } else if (numSeparators == 2) {


    // First, place all non VE separators first in the order, since
    // these will be affected by the pruning parameters (while the
    // VE seps (if any) don't change with pruning.
    if (section.separators[clique.ceReceiveSeparators[0]].veSeparator
	&& !section.separators[clique.ceReceiveSeparators[1]].veSeparator) {
      // sep 0 is a VE separator, but sep 1 isn't. Put the VE sep last.
      swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[1]);
    } else if (section.separators[clique.ceReceiveSeparators[1]].veSeparator
	       && !section.separators[clique.ceReceiveSeparators[0]].veSeparator) {
      // sep 1 is a VE separator, but sep 0 isn't, do nothing.
    } else {
      // If neither or both of the separators are VE separators, then
      // iterate through smaller weight separator first. If no
      // intersection, then this still can have a benefit since we
      // want to do the most work in the inner (rather than the outer)
      // loops. Of course, we don't know the real cost of the
      // separater state space until run time (which depends on the
      // evidence set, current parameters, and the pruning
      // thresholds), but we use the weight heuristic for now.
      if (section.separators[clique.ceReceiveSeparators[0]].weight() <=
	  section.separators[clique.ceReceiveSeparators[1]].weight()) {
	// do nothing
      } else {
	swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[1]);
      }
    }


    // shortcuts to separator 0 and 1
    SeparatorClique& s0 = section.separators[clique.ceReceiveSeparators[0]];
    SeparatorClique& s1 = section.separators[clique.ceReceiveSeparators[1]];


    // intersection of the two separators
    set<RV*> sepIntersection;
    set_intersection(s0.nodes.begin(),s0.nodes.end(),
		     s1.nodes.begin(),s1.nodes.end(),
		     inserter(sepIntersection,sepIntersection.end()));

    // first one in order should be empty.
    s0.accumulatedIntersection.clear();
    // and remainder of first one should get the rest
    s0.remainder = s0.nodes;
    // 2nd one in should be intersection
    s1.accumulatedIntersection = sepIntersection;
    // and remainder gets the residual.
    set_difference(s1.nodes.begin(),
		   s1.nodes.end(),
		   sepIntersection.begin(),sepIntersection.end(),
		   inserter(s1.remainder,
			    s1.remainder.end()));
    // compute union of all separators.
    set_union(s0.nodes.begin(),s0.nodes.end(),
	      s1.nodes.begin(),s1.nodes.end(),
	      inserter(clique.unionIncommingCESeps,clique.unionIncommingCESeps.end()));

    assert ( s0.accumulatedIntersection.size() + s0.remainder.size() == s0.nodes.size() );
    assert ( s1.accumulatedIntersection.size() + s1.remainder.size() == s1.nodes.size() );

  } else {
    // there are 3 or more separators, determine proper order and then
    // compute running accumulated intersection relative to that order.

    // Goal: to produce an optimal separator iteration order.  In
    // general, the number of "live" surviving entries of all
    // separators when they are intersected the end of a sep
    // traversal will be the same no matter the order that is
    // chosen. What we want here, however, is something that kills
    // non-live entries ASAP, to avoid extra work.  This will depend
    // on properties of the separators (which we have here) and the
    // entries contained within the separators (which we don't have
    // until run time). Therefore, our order is a heuristic.
    // We do this by reordering the entries in clique.ceReceiveSeparators.

    // Other Ideas:
    // Probably a dynamic programing algorithm would work better
    // here. Use DP to find the order that minimizes the sum (in
    // reverse order) of the intersection of the last set with all
    // previous ones.

    //
    // What we first do: place VE seps last (since they are uneffected by pruning)
    // What we next do: sort in reverse order, by separator that has the
    // minimum intersection with others. At each step, when we have
    // found the sep with min intersection among current number of
    // seps, we place this at the *end* of the separation iteration
    // order. Once we get back down to two separators, we start with
    // the one of minimal weight (to do as much as possible in the
    // inner loops).  Also, when ties exist, do the thing that
    // causes the fewest number of branches and has the inner most
    // loop run interuppted for the greatest amount of time (i.e.,
    // min weight goes first).

    // printf("before sorting\n");
    // for (unsigned i = 0; i< clique.ceReceiveSeparators.size(); i++) {
    // printf("ceRecSep[%d] = %d\n",i,clique.ceReceiveSeparators[i]);
    // }

    // TODO: need to change this to be based only on hidden nodes,
    // since the intersection of observed nodes doesn't count
    // (observed nodes don't contribute to intersection pruning)
    // since they are never inserted into the clique.

#if 0
    for (unsigned i=0;i<clique.ceReceiveSeparators.size();i++) {
      if ((clique.ceReceiveSeparators[i] == (section.separators.size()-1))) {
	swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[i]);
	break;
      }
    }
#endif

    // First, place all non VE separators first in the order, since
    // these will be affected by the pruning parameters (while the
    // VE seps (if any) don't change with pruning.
    int swapPosition = numSeparators-1;
    for (int i=0;i<swapPosition;) {
      if (section.separators[clique.ceReceiveSeparators[i]].veSeparator) {
	swap(clique.ceReceiveSeparators[i],clique.ceReceiveSeparators[swapPosition]);
	swapPosition--;
      } else
	i++;
    }
    if (!section.separators[clique.ceReceiveSeparators[swapPosition]].veSeparator)
      swapPosition++;
    const int firstVESeparator = swapPosition;
    const int lastRealSeparator = swapPosition-1;

    // printf("lastReal = %d, ceseps:",lastRealSeparator); 
    // for (unsigned i=0;i<numSeparators;i++) {
    // printf("%d ",clique.ceReceiveSeparators[i]);
    // }
    // printf("\n");

    // printf("Sorting Seps\n");

    // next, sort the real separators based on maximizing intersection.
    if (lastRealSeparator > 0) {
      // then we have at least 2 real separators.

      unsigned lastSeparator = lastRealSeparator;
      count_iterator<set <RV*> > sep_intr_set;
      set<RV*> sep_union_set;
      set<RV*> empty;
      vector < pair<unsigned,unsigned> > sepIntersections; 
      while (lastSeparator > 1) {
	// then at least three real separators to go.

	// compute separator index with minimum intersection with all previous separators
	// and swap its index with lastSeparator.
	sepIntersections.clear();
	sepIntersections.resize(lastSeparator+1);
	for (unsigned i=0;i<=lastSeparator;i++) {
	  SeparatorClique& sep_i =
	    section.separators[clique.ceReceiveSeparators[i]];
	  sepIntersections[i].second = i;
	  sepIntersections[i].first = 0;
	  sep_union_set.clear();
	  for (unsigned j=0;j<=lastSeparator;j++) {
	    if (j == i)
	      continue;
	    SeparatorClique& sep_j =
	      section.separators[clique.ceReceiveSeparators[j]];
	    set_union(empty.begin(),empty.end(),
		      sep_j.nodes.begin(),sep_j.nodes.end(),
		      inserter(sep_union_set,sep_union_set.end()));
	  }
	  sep_intr_set.reset();
	  set_intersection(sep_i.nodes.begin(),sep_i.nodes.end(),
			   sep_union_set.begin(),sep_union_set.end(),
			   sep_intr_set);
	  sepIntersections[i].first = sep_intr_set.count();
	}

	// sort in (default) ascending order
	sort(sepIntersections.begin(),sepIntersections.end(),PairUnsigned1stElementCompare());

	// among number of ties (the ones with minimal intersection
	// with others), find the one with greatest weight to place
	// last.
	unsigned bestIndex = 0;
	float bestWeight = section.separators[clique.ceReceiveSeparators[sepIntersections[bestIndex].second]].weight();
	unsigned i = 1;
	while (i < sepIntersections.size() && sepIntersections[bestIndex].first == sepIntersections[i].first) {
	  float curWeight = section.separators[clique.ceReceiveSeparators[sepIntersections[i].second]].weight();
	  if (curWeight > bestWeight) {
	    bestWeight = curWeight;
	    bestIndex = i;
	  }
	  i++;
	}

	// finally, swap the best into last place.
	// printf("swaping sep %d with %d\n",bestIndex,lastSeparator);
	swap(clique.ceReceiveSeparators[sepIntersections[bestIndex].second],clique.ceReceiveSeparators[lastSeparator]);
	lastSeparator --;
      }

      // Lastly, the first two separators should be iterated
      // so that the min weight one goes first.
      if (section.separators[clique.ceReceiveSeparators[0]].weight() >
	  section.separators[clique.ceReceiveSeparators[1]].weight()) {
	swap(clique.ceReceiveSeparators[0],clique.ceReceiveSeparators[1]);
      }
    }

    // next, sort  the VE separators to maximize overlap.
    const unsigned numVESeparators = numSeparators-firstVESeparator;
    if (numVESeparators > 1) {
      // so at least two VE separators.

      int lastSeparator = numSeparators-1;
      count_iterator<set <RV*> > sep_intr_set;
      set<RV*> sep_union_set;
      set<RV*> empty;
      vector < pair<unsigned,unsigned> > sepIntersections; 
      while (lastSeparator > firstVESeparator + 1) {
	// then at least three to go.

	// compute separator index with minimum intersection with
	// *all* previous separators (both normal and VE) and swap
	// its index with lastSeparator.
	sepIntersections.clear();
	sepIntersections.resize(lastSeparator+1-firstVESeparator);
	for (int i=firstVESeparator;i<=lastSeparator;i++) {
	  SeparatorClique& sep_i =
	    section.separators[clique.ceReceiveSeparators[i]];
	  sepIntersections[i-firstVESeparator].second = i;
	  sepIntersections[i-firstVESeparator].first = 0;
	  sep_union_set.clear();
	  for (int j=0;j<=lastSeparator;j++) {
	    if (j == i)
	      continue;
	    SeparatorClique& sep_j =
	      section.separators[clique.ceReceiveSeparators[j]];
	    set_union(empty.begin(),empty.end(),
		      sep_j.nodes.begin(),sep_j.nodes.end(),
		      inserter(sep_union_set,sep_union_set.end()));

	  }
	  sep_intr_set.reset();
	  set_intersection(sep_i.nodes.begin(),sep_i.nodes.end(),
			   sep_union_set.begin(),sep_union_set.end(),
			   sep_intr_set);
	  sepIntersections[i-firstVESeparator].first = sep_intr_set.count();
	}


	// 	  printf("before sorting: VE seps have intersection among others:");
	// 	  for (int i=firstVESeparator;i<=lastSeparator;i++) {
	// 	    printf("%d has %d,",
	// 		   sepIntersections[i-firstVESeparator].second,
	// 		   sepIntersections[i-firstVESeparator].first);
	// 	  }
	// 	  printf("\n");


	// sort in (default) ascending order
	sort(sepIntersections.begin(),sepIntersections.end(),PairUnsigned1stElementCompare());

	//  	  printf("after sorting: VE seps have intersection among others size = %d:",sepIntersections.size());
	//  	  for (int i=firstVESeparator;i<=lastSeparator;i++) {
	//  	    printf("Sep %d w inter %d.  ",
	//  		   clique.ceReceiveSeparators[sepIntersections[i-firstVESeparator].second],
	//  		   sepIntersections[i-firstVESeparator].first);
	//  	  }
	// 	  printf("\n");

	// 	  vector < pair<unsigned,unsigned> >::iterator it;
	// 	  for (it = sepIntersections.begin(); 
	// 		 ((it) != sepIntersections.end()); it++) {
	// 	    printf("cur=(%d,%d), next=(%d,%d), cur<next=%d\n",
	// 		   (*it).first,(*it).second,
	// 		   (*(it+1)).first,(*(it+1)).second,
	// 		   (*it)<(*(it+1)));
	// 	  }


	// among number of ties (the ones with minimal intersection
	// with others), find the one with greatest weight to place
	// last.
	unsigned bestIndex = 0;
	float bestWeight = section.separators[clique.ceReceiveSeparators[sepIntersections[bestIndex].second]].weight();
	unsigned i = 1;
	while (i < sepIntersections.size() && sepIntersections[bestIndex].first == sepIntersections[i].first) {
	  float curWeight = section.separators[clique.ceReceiveSeparators[sepIntersections[i].second]].weight();
	  if (curWeight > bestWeight) {
	    bestWeight = curWeight;
	    bestIndex = i;
	  }
	  i++;
	}

	// finally, swap the best into last place.
	// printf("swaping sep %d with %d\n",bestIndex,lastSeparator);
	swap(clique.ceReceiveSeparators[sepIntersections[bestIndex].second],clique.ceReceiveSeparators[lastSeparator]);
	lastSeparator --;
      }

      // Lastly, the first two separators should be iterated
      // so that the min weight one goes first.
      if (section.separators[clique.ceReceiveSeparators[firstVESeparator]].weight() >
	  section.separators[clique.ceReceiveSeparators[firstVESeparator+1]].weight()) {
	swap(clique.ceReceiveSeparators[firstVESeparator],clique.ceReceiveSeparators[firstVESeparator+1]);
      }
    }


    ///////////////////////////////////////////////////////
    // DONE SORTING clique.ceReceiveSeparators
    // printf("after sorting\n");
    // for (unsigned i = 0; i< clique.ceReceiveSeparators.size(); i++) {
    // printf("ceRecSep[%d] = %d\n",i,clique.ceReceiveSeparators[i]);
    // }
    // Compute the cummulative intersection of the sepsets using the
    // currently assigned sepset order.
    {

      // initialize union of all previous separators

      clique.unionIncommingCESeps = 
	section.separators[clique.ceReceiveSeparators[0]].nodes;
      section.separators[clique.ceReceiveSeparators[0]].accumulatedIntersection.clear();
      section.separators[clique.ceReceiveSeparators[0]].remainder
	= section.separators[clique.ceReceiveSeparators[0]].nodes;

      for (unsigned sep=1;sep<numSeparators;sep++) {

	// reference variables for easy access
	set<RV*>& sepNodes
	  = section.separators[clique.ceReceiveSeparators[sep]].nodes;
	set<RV*>& sepAccumInter
	  = section.separators[clique.ceReceiveSeparators[sep]].accumulatedIntersection;


	// create the intersection of 1) the union of all previous nodes in
	// the sep order, and 2) the current sep nodes.
	sepAccumInter.clear();
	set_intersection(clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
			 sepNodes.begin(),sepNodes.end(),
			 inserter(sepAccumInter,sepAccumInter.end()));

	// compute the separator remainder while we're at it.
	// specifically: remainder = sepNodes - sepAccumInter.
	section.separators[clique.ceReceiveSeparators[sep]].remainder.clear();
	set_difference(sepNodes.begin(),sepNodes.end(),
		       sepAccumInter.begin(),sepAccumInter.end(),
		       inserter(section.separators[clique.ceReceiveSeparators[sep]].remainder,
				section.separators[clique.ceReceiveSeparators[sep]].remainder.end()));


	// update the accumulated (union) of all previous sep nodes.
	set<RV*> res;
	set_union(sepNodes.begin(),sepNodes.end(),
		  clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
		  inserter(res,res.end()));
	clique.unionIncommingCESeps = res;	
      }
    }

  }


  // lastly, assign unassignedIteratedNodes in this clique
  {
    // unassigned nodes, unassigned iterated nodes, 
    // compute: unassignedIteratedNodes  = nodes - (assignedNodes U unionIncommingCESeps)
    // first: compute res = nodes - assignedNodes
    set<RV*> res;
    set_difference(clique.nodes.begin(),clique.nodes.end(),
		   clique.assignedNodes.begin(),clique.assignedNodes.end(),
		   inserter(res,res.end()));
    // next: compute unassignedIteratedNodes = res - unionIncommingCESeps
    // note at this point unionIncommingCESeps contains the union of all nodes in all separators
    set_difference(res.begin(),res.end(),
		   clique.unionIncommingCESeps.begin(),clique.unionIncommingCESeps.end(),
		   inserter(clique.unassignedIteratedNodes,
			    clique.unassignedIteratedNodes.end()));
  }
  // @@@ perhaps subtract out unassignedInSection since they are different??
}



/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::computeSeparatorIterationOrders()
 *   computes the order in which we iterate over the separators
 *   when we do seperator driven clique potential function
 *   creation.
 *
 * Preconditions:
 *   separators both within and between sections must have been
 *   created (meaning createSeparators() must be called).
 *
 * Postconditions:
 *   Separator orders have now been created.
 *
 *
 * Side Effects:
 *    Changes member variables within cliques of each section.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::computeSeparatorIterationOrders() {
  computeSeparatorIterationOrders(P1);
  computeSeparatorIterationOrders(Co);
  computeSeparatorIterationOrders(E1);
}



// compute the assignment order for nodes in this
// section's cliques relative to each clique's incomming separators, and while
// doing so, also set the dispositions for each of the resulting
// nodes in each clique.


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions()
 *
 *   Computes for each clique the union of the set of nodes which have
 *   been iterated unassigned in previous cliques in the JT. This is
 *   used to determine which, in each clique, nodes should be iterated
 *   over and which are already assigned by a separator driven
 *   iteration.
 *
 * Preconditions:
 *   computeSeparatorIterationOrder() must have been called. Section
 *   in argument 'section' must not be empty.
 *
 * Postconditions:
 *   cliques know which nodes to iterate over or not.
 *
 * Side Effects:
 *    Changes member variables within cliques of this section.
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& section, 
								char const *varCliqueAssignmentPrior)
{
  if (varCliqueAssignmentPrior && strlen(varCliqueAssignmentPrior) > 0) {
    infoMsg(IM::Inference, IM::Giga,"Sorting cliques variables using priority order (%s)\n",varCliqueAssignmentPrior);
  }
  for (unsigned i=0; i < section.cliques.size(); i++) {
    section.cliques[i].sortAndAssignDispositions(varCliqueAssignmentPrior);
    if (IM::messageGlb(IM::Inference, IM::Max)) {
      printf("Clique %d variables after sort:",i);
      printRVSet(stdout,section.cliques[i].sortedAssignedNodes,true);
    }
  }
}

void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions(char const *varCliqueAssignmentPrior) {
  // printf("sorting P1\n");
  sortCliqueAssignedNodesAndComputeDispositions(P1,varCliqueAssignmentPrior);
  // printf("sorting Co\n");
  sortCliqueAssignedNodesAndComputeDispositions(Co,varCliqueAssignmentPrior);
  // printf("sorting E1\n");
  sortCliqueAssignedNodesAndComputeDispositions(E1,varCliqueAssignmentPrior);
}


void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions(JT_Partition& section) {
  for (unsigned i=0; i < section.cliques.size(); i++) {
    section.cliques[i].sortAndAssignDispositions();
    if (IM::messageGlb(IM::Inference, IM::Max)) {
      printf("Clique %d variables after sort:",i);
      printRVSet(stdout,section.cliques[i].sortedAssignedNodes,true);
    }
  }
}

void 
SectionScheduler::sortCliqueAssignedNodesAndComputeDispositions() {
  // printf("sorting P1\n");
  sortCliqueAssignedNodesAndComputeDispositions(P1);
  // printf("sorting Co\n");
  sortCliqueAssignedNodesAndComputeDispositions(Co);
  // printf("sorting E1\n");
  sortCliqueAssignedNodesAndComputeDispositions(E1);
}

 
/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::setUpJTDataStructures()
 *
 *   Sets up all the data structures in a JT (other than preparing for
 *   unrolling), and calls all the routines in the necessary
 *   order. This is like the main routine of this class, as it calls
 *   what needs to be done to set up a JT.
 *
 * Preconditions:
 *   Junction tree must have been created. Should not have called 
 *   setUpDataStructures() before. The C sections in the template
 *   have variables, but the other sections (P, and E) might be empty.
 *
 * Postconditions:
 *   All data structures are set up. The JT is ready to have
 *   prepareForUnrolling() called.
 *
 * Side Effects:
 *   Changes many internal data structures in this object. 
 *
 * Results:
 *    None.
 *
 *-----------------------------------------------------------------------
 */
void
SectionScheduler::setUpJTDataStructures(const char* varSectionAssignmentPrior,
					const char *varCliqueAssignmentPrior)
{
  // main() routine for this class.

  createSectionJunctionForests();   // create undirected graph with edges in clique[i].neighbors

  computeSectionInterfaces();       // setup left & right interfaces
  //createFactorCliques();
  createDirectedGraphOfCliques();   // create junction forest with edges in clique[i].children
  //  assignRVsToCliques(varSectionAssignmentPrior,varCliqueAssignmentPrior);
  //assignFactorsToCliques();
  //  computeUnassignedCliqueNodes();
  // TODO: assignScoringFactorsToCliques();
  //setUpMessagePassingOrders();
  // create seps and VE seps.
  createSeparators();               // create within- & between-section separators per message order

  computeSeparatorIterationOrders();

  // -- -- used only to compute weight.
  getCumulativeUnassignedIteratedNodes(); 
  // -- --

sortCliqueAssignedNodesAndComputeDispositions();
  //sortCliqueAssignedNodesAndComputeDispositions(varCliqueAssignmentPrior);
}





/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::junctionTreeWeight()
 *
 * Compute the 'junction tree weight' (roughly, the log10(cost of
 * doing inference)) for the set of cliques given in cliques. Note,
 * cliques *must* be a valid set of maxcliques of a junction tree --
 * if they are not, unexpected results are returned.
 *
 * Note also, this returns an upper bound on the cost of a JT in this
 * section, rather than the actual cost. It returns an upper bound
 * since it assumes that the separator driven iteration will be at
 * full cardinality (it can't tell the difference between previous
 * assigned nodes which have been heavily cut back and previous
 * unassigned nodes which are really at full cardinality). Therefore,
 * this routine returns the conservative upper bound. The goal
 * is to minimize this upper bound.
 *
 *
 * Preconditions:
 *   The section must have been fully instantiated. I.e.,
 *   we must have that assignRVsToCliques have been called.
 *
 * Postconditions:
 *   The weight, as inference would do it, of this clique is 
 *   computed.
 *
 * Side Effects:
 *   none.
 *
 * Results:
 *   the weight
 *
 *-----------------------------------------------------------------------
 */
double
SectionScheduler::junctionTreeWeight(JT_Partition& part,
				     const unsigned rootClique,
				     vector< set<RV*> > *lp_nodes,
				     vector< set<RV*> > *rp_nodes)
{
  // check for case of empty E or P section.
  if (part.cliques.size() == 0)
    return 0;

  // FIXME: make this the sum of the now multiple root-ish cliques?

  MaxClique& curClique = part.cliques[rootClique];

  // @@ add in unassigned in section information to next call
  double weight = curClique.weightInJunctionTree(part.unassignedInPartition,
						 jtWeightUpperBound,
						 jtWeightMoreConservative,
						 true,
						 lp_nodes,rp_nodes);
  // double weight = curClique.weight();
  for (unsigned childNo=0;
       childNo<curClique.children.size();childNo++) {
    double child_weight = junctionTreeWeight(part,
					     curClique.children[childNo],lp_nodes,rp_nodes);
    // i.e., log addition for weight = weight + child_weight
    weight = log10add(weight,child_weight);
  }
  return weight;
}


/*-
 *-----------------------------------------------------------------------
 * SectionScheduler::junctionTreeWeight()
 *
 * Given a set of maxcliques for a section, and an interface for
 * this (can be left right, or any set including empty, the only
 * condition is that it must be covered by at least one of the
 * cliques), compute the junction tree for this set and return the
 * estimated JT cost. This is a static routine so can be called from
 * anywhere.
 *
 * Preconditions:
 *   one of the cliques must cover interface nodes. I.e.,
 *   there must be an i such that cliques[i].nodes >= interfaceNodes.
 *
 * Postconditions:
 *   An estimate of the weight, as in inference, of this clique set is 
 *   computed and returned.
 *
 * Side Effects:
 *   none.
 *
 * Results:
 *   the weight
 *
 *-----------------------------------------------------------------------
 */
double
SectionScheduler::junctionTreeWeight(vector<MaxClique>& cliques,
				     const vector < set<RV*> > &interfaceNodes,
				     vector< set<RV*> > *lp_nodes,
				     vector< set<RV*> > *rp_nodes)

{

  // const bool is_P_section = (lp_nodes == NULL);
  const bool is_E_section = (rp_nodes == NULL);

  // might be a section that is empty.
  if (cliques.size() == 0)
    return 0.0;

  Section part;
  const set <RV*> emptySet;
  const vector< set <RV*> > emptyInterface;
  part.cliques = cliques;


  // set up the nodes into part
  for (unsigned i=0;i<cliques.size();i++) {
    // while we're at it, clear up the JT structures we're
    // about to create.
    cliques[i].clearJTStructures();
    set_union(emptySet.begin(),emptySet.end(),
	      cliques[i].nodes.begin(),cliques[i].nodes.end(),
	      inserter(part.nodes,part.nodes.end()));
  }
  createSectionJunctionTree(part);

  // a JT version of this section.
  JT_Partition jt_part(part,emptyInterface,interfaceNodes);

  bool tmp;
  vector<unsigned> root( max(1U, (unsigned)interfaceNodes.size()) );
  if (!is_E_section) 
    jt_part.findRInterfaceClique(root,tmp,interfaceCliquePriorityStr);
  else {
    // 'E_root_clique case'
    // Presumably, this is an E section, so the root should be done
    // same as E_root_clique computed above.
    // "update E_root_clique"
    // root = 0; 
    // currently, find max weight clique.
    double weight = -10e20;
    for (unsigned i=0;i<cliques.size();i++) {
      if (MaxClique::computeWeight(cliques[i].nodes) > weight) {
	// FIXME - should this be the sum of the max weights in each connected component?
	root[0] = i;
      }
    }
  }

  createDirectedGraphOfCliques(jt_part,root[0]);
  assignRVsToCliques("candidate section",jt_part,root[0],"","");

  vector< pair<unsigned,unsigned> > message_order;
  vector< unsigned > leaf_cliques;
  setUpMessagePassingOrder(jt_part,
			   root[0],
			   message_order,
			   ~0x0,
			   leaf_cliques);

  // Note that since we're calling createSeparators(jt_part,msg_order)
  // here directly, this means that the left interface clique will not
  // have its separator to the left section filled in.
  createSeparators(jt_part,message_order);

  computeSeparatorIterationOrders(jt_part);
  getCumulativeUnassignedIteratedNodes(jt_part,root[0]);

  // If the code is changed so that jt weight needs the sorted
  // assigned nodes or the dispositions, we'll need to make a static version of
  // the following routine and call it here.
  // sortCliqueAssignedNodesAndComputeDispositions();

  // return jt_part.cliques[root].cumulativeUnassignedIteratedNodes.size();
  double weight = junctionTreeWeight(jt_part,root[0],lp_nodes,rp_nodes);

  if (jtWeightPenalizeUnassignedIterated > 0.0) {
    unsigned badness_count=0;
    set <RV*>::iterator it;
    for (it = jt_part.cliques[root[0]].cumulativeUnassignedIteratedNodes.begin();
	 it != jt_part.cliques[root[0]].cumulativeUnassignedIteratedNodes.end();
	 it++) 
      {
	RV* rv = (*it);

	assert ( rv-> discrete() );
	DiscRV* drv = (DiscRV*)rv;
	if (!drv->sparse())
	  continue;

	// if it has not at all been assigned in this section, we
	// assume it has been assigned in previous section, and dont'
	// count it's cost.
	if (jt_part.unassignedInPartition.find(rv) ==
	    jt_part.unassignedInPartition.end())
	  continue;

	badness_count ++;
      }
    weight = badness_count*jtWeightPenalizeUnassignedIterated + weight;
  }
  return weight;

  // return badness_count*1000 + weight;
  // return badness_count;

}

