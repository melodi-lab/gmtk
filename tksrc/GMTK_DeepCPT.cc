/*-
 *
 * GMTK_DeepCPT.cc
 *     Deep neural network CPT.
 *
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes 
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *   This CPT gives GMTK the ability to use a deep neural netowrk
 *   to compute probabilities.
 *
 *   This CPT does P(X | \Pi_X = j) = NN(I(j)) where
 *     - \Pi_X are the parents of RV X
 *     - NN is a DeepNN with |Domain(X)| outputs and \Sum_{Y \in \Pi_X} |Domain(Y)| inputs
 *     - I produces an |\Pi_X|-hot vector (0s everywhere except 1s in the positions
 *       corresponding to the value(s) observed for \Pi_X) of length \Sum_{Y \in \Pi_X} |Domain(Y)| 
 *
 * An example of how you'd specify DeepCPTs:
 * (in a masterfile:)


DEEP_CPT_IN_FILE inline 1
0
cptName
3             % number of parents
4 2 6         % parent cardinalities, total must match DeepNN's # inputs
2             % self cardinality 
deepNNName    % DeepNN to use to compute probabilities, # outputs must
              %   match self cardinality


 *
 *
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <typeinfo>

#include "general.h"
#include "error.h"
#include "rand.h"
#include "phipac.h"

#include "GMTK_DeepCPT.h"
#include "GMTK_GMParms.h"
#include "GMTK_DiscRV.h"

VCID(HGID)


////////////////////////////////////////////////////////////////////
//        General create, read, destroy routines 
////////////////////////////////////////////////////////////////////



/*-
 *-----------------------------------------------------------------------
 * DeepCPT::read(is)
 *      read in the table.
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Configures the DeepCPT instance
 *
 *-----------------------------------------------------------------------
 */


void
DeepCPT::read(iDataStreamFile& is)
{
  string str;

  ///////////////////////////////////////////////////////////////////////
  // The first set of options match that of any other CPT. Namely,
  // we have
  //   1) a name, 
  //   2) num parents
  //   3) parent cardinalities
  //   4) self cardinaltiy
  //   5) DeepNN name to use to get probabilities

  // read the name of the object.
  NamedObject::read(is);
  is.read(_numParents,"Can't read DeepCPT's num parents");
  if (_numParents < 1) 
    error("ERROR: reading file '%s' line %d, DeepCPT '%s' trying to use %d < 1 num parents.",
	  is.fileName(),is.lineNo(),name().c_str(),_numParents);
  if (_numParents >= warningNumParents)
    warning("WARNING: creating DeepCPT with %d parents in file '%s' line %d",
	    _numParents,is.fileName(),is.lineNo());
  cardinalities.resize(_numParents);

  // read the parent cardinalities
  unsigned num_inputs = 0;
  for (unsigned i=0;i<_numParents;i++) {
    is.read(cardinalities[i],"Can't read DeepCPT parent cardinality");
    if (cardinalities[i] <= 0)
      error("ERROR: reading file '%s' line %d, DeepCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	    is.fileName(),is.lineNo(),name().c_str(),cardinalities[i],i);
    num_inputs += cardinalities[i];
  }
  curParentValues.resize(_numParents);
  cachedParentValues.resize(_numParents);

  // read the self cardinality
  is.read(_card,"Can't read DeepCPT self cardinality");
  if (_card <= 0)
    error("ERROR: reading file '%s' line %d, DeepCPT '%s' trying to use 0 or negative (%d) cardinality table, position %d.",
	  is.fileName(),is.lineNo(),name().c_str(),_card,_numParents);

  cached_CPT = NULL;

  is.read(str, "Can't read DeepCPT's DeepNN name");
  if (GM_Parms.deepNNsMap.find(str) == GM_Parms.deepNNsMap.end()) {
    error("ERROR: reading file '%s' line %d,  DeepCPT '%s' uses undefined DeepNN '%s'\n", is.fileName(), is.lineNo(), name().c_str(), str.c_str());
  }
  dmlp = GM_Parms.deepNNs[ GM_Parms.deepNNsMap[str] ];
  assert(dmlp);

  // still here? Do more error checking.
  
  // check that layer 0 matches the input vector size (+1 for bias element)
  if (dmlp->numInputs() != num_inputs) {
    error("ERROR: reading file '%s' line %d, DeepCPT '%s' uses DeepNN with %u inputs, but it's parent cardinalities sum to %d\n", is.fileName(), is.lineNo(), name().c_str(), dmlp->numInputs(), num_inputs);
  }

  // check that the final # outputs = parent cardinality
  if (dmlp->numOutputs() != _card) {
    error("ERROR: reading file '%s' line %d, DeepCPT '%s' requires DeepNN '%s' to have %u outputs, but it has %u\n",
	  is.fileName(), is.lineNo(), _name.c_str(), dmlp->name().c_str(), _card, dmlp->numOutputs());
  }
  setBasicAllocatedBit();
}



/*-
 *-----------------------------------------------------------------------
 * DeepCPT::write(os)
 *      write out data to file 'os'. 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      No effects other than  moving the file pointer of os.
 *
 *-----------------------------------------------------------------------
 */
void
DeepCPT::write(oDataStreamFile& os)
{
  assert ( basicAllocatedBitIsSet() );
  NamedObject::write(os);
  os.nl();
  os.write(_numParents,"DeepCPT::write numParents"); 
  os.writeComment("number parents");os.nl();
  for (unsigned i=0;i<_numParents;i++) {
    os.write(cardinalities[i],"DeepCPT::write cardinality");
  }
  os.write(card(),"DeepCPT::write cardinality");
  os.writeComment("cardinalities");
  os.write(dmlp->name()); os.nl();
  os.writeComment("DeepNN name");
  os.nl();
}


////////////////////////////////////////////////////////////////////
//        Probability routines (see GMTK_CPT.h for docs)
////////////////////////////////////////////////////////////////////


void 
DeepCPT::becomeAwareOfParentValues( vector <RV *>& parents, const RV* rv ) {
  assert ( basicAllocatedBitIsSet() );
  assert ( parents.size() == _numParents );
  
  for (unsigned i = 0; i < _numParents; i++) {
    if (((DiscRV*)parents[i])->val >= cardinalities[i])
      error("ERROR: DeepCPT %s used by RV %s(%d), invalid parent value for parent %s(%d) (parent number %d), parentValue = %d but parent RV cardinality = %d\n",
	    name().c_str(),
	    rv->name().c_str(),rv->frame(),
	    parents[i]->name().c_str(),parents[i]->frame(),
	    i,
	    ((DiscRV*)parents[i])->val,cardinalities[i]);

    curParentValues[i] = RV2DRV(parents[0])->val;
  }
}


logpr
DeepCPT::applyNN(DiscRV * drv) {
  register DiscRVType val = drv->val;

  logpr p((void*)NULL);

  // check if the CPT is cached
  bool hit = cached_CPT != NULL;
  for (unsigned i=0; i < _numParents && hit; i+=1)
    hit = hit && (curParentValues[i] == cachedParentValues[i]);
  if (hit) {
    p.setFromP(cached_CPT[val]);
    return p;
  }

  // Not in the cache, so compute & cache it

  cachedParentValues = curParentValues;

  // assemble input vector

  unsigned num_inputs = dmlp->numInputs() + 1;
  float *input_vector = new float[num_inputs];
  memset(input_vector, 0, num_inputs * sizeof(float));

  input_vector[num_inputs-1] = 1.0; // homogeneous coordinates
  unsigned offset = 0;
  for (unsigned i=0; i < _numParents; i+=1) {
    input_vector[offset + curParentValues[i]] = 1.0;
    offset += cardinalities[i];
  }

  if (cached_CPT) delete[] cached_CPT;
  cached_CPT = dmlp->applyDeepModel(input_vector);

  // log the CPT entry
  p.setFromP(cached_CPT[val]);
  return p;
}

logpr 
DeepCPT::probGivenParents(vector <RV *>& parents, DiscRV * drv) {
  assert ( bitmask & bm_basicAllocated );
  DeepCPT::becomeAwareOfParentValues(parents, drv);
  return applyNN(drv);
}


void DeepCPT::becomeAwareOfParentValuesAndIterBegin(vector < RV *>& parents, 
						  iterator &it,
						  DiscRV* drv,
						  logpr& p) 
{
  DeepCPT::becomeAwareOfParentValues(parents,drv);
  begin(it,drv,p);
}

void DeepCPT::begin(iterator& it,DiscRV* drv, logpr& p) {
  assert ( bitmask & bm_basicAllocated );
  it.setCPT(this);
  it.drv = drv;
  drv->val = 0;
  do {
    p = applyNN(drv);
    if (! p.essentially_zero()) break;
    drv->val++;
    // We keep the following check as we must have that at least one
    // entry is non-zero.  The read code of the MDCPT should ensure
    // this as sure all parameter update procedures, as long as
    // normalizationThreshold is not set to large.
    if (drv->val >= card()) {
      error("ERROR: DeepCPT '%s' used for RV '%s(%d)' has a DeepNN outout with all zeros. Program Exiting ...\n",
	    name().c_str(),drv->name().c_str(),drv->frame());
    }
  } while (1);
}

bool DeepCPT::next(iterator &it,logpr& p) {
  do {
    it.drv->val++;
    // don't increment past the last value.
    if (it.drv->val == card()) {
      return false;
    }
    p = applyNN(it.drv);
  } while (p.essentially_zero());
  return true;
}



/*-
 *-----------------------------------------------------------------------
 * DeepCPT::random sample(is)
 *      random sample according to current 'distribution' 
 * 
 * Results:
 *      No results.
 *
 * Side Effects:
 *      Changes val
 *
 *-----------------------------------------------------------------------
 */
int DeepCPT::randomSample(DiscRV* drv)
{
  // sum up the observations in current frame.
  // TODO: this routine is currently not correct (but it is not used).
#if 0
  logpr sum;
  for (unsigned i=0;i<obs->numContinuous();i++) {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
    sum += tmp;
  }
  logpr uniform = (rnd.drand48()*sum.unlog());
  sum.set_to_zero();
  unsigned i=0;
  do {
    logpr tmp((void*)NULL);
    tmp.valref() = (*obs->floatVecAtFrame(drv->frame(),i));
    sum += tmp;
    if (uniform <= sum)
      break;
  } while (++i < obs->numContinuous());

  if (veMode == VE_Dense) {
    return i;
  } else {
    // this assumes that ints in sparse case are ordered
    return *(obs->unsignedVecAtFrame(drv->frame()) + i); 
  }
#endif
  return 0;
}



////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN

#endif
