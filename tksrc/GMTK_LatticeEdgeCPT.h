/*-
 * GMTK_LatticeEdgeCPT.h
 *           CPT for lattice edge misstions
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


#ifndef GMTK_LATTICE_EDGE_CPT_H
#define GMTK_LATTICE_EDGE_CPT_H


#include "GMTK_CPT.h"
#include "GMTK_LatticeADT.h"


/**
 * lattice edge CPT
 */
class LatticeEdgeCPT : public CPT {
public:
	LatticeEdgeCPT();
	~LatticeEdgeCPT();

	//////////////////////////////////////////////////////////
	// virtual functions from class CPT
                                                                                
        // Allocate memory, etc. for the internal data structures
        // for this CPT, depending on current _numParents & cardinalities.
        virtual void allocateBasicInternalStructures() {}
                                                                                
        ///////////////////////////////////////////////////////////
        // Probability evaluation, compute Pr( child | parents ), and
        // iterator support. See GMTK_CPT.h for documentation.
        virtual void becomeAwareOfParentValues(vector< RV* >& parents, const RV* rv) { assert(0); }
        virtual void begin(iterator& it, DiscRV* drv,logpr& p) { assert(0); }
        virtual void becomeAwareOfParentValuesAndIterBegin(vector< RV* >& parents, iterator &it, DiscRV* drv, logpr& p);
        virtual logpr probGivenParents(vector< RV* >& parents, DiscRV* drv);
        virtual bool next(iterator &it,logpr& p);
 
        // Include here an extra routine that returns the probability
        // of the child 'val' given the parents are the assigned to
        // the set of values set during the most previous call to
        // becomeAwareOfParentValues.
        logpr probGivenParents(const int val) { assert(0); }
 
        ///////////////////////////////////////////////////////////
        // Given the current parent values, generate a random sample.
        virtual int randomSample(DiscRV*drv) { assert(0); }
 
        ///////////////////////////////////////////////////////////
        // Re-normalize the output distributions
        virtual void normalize() {}
        // set all values to random values.
        virtual void makeRandom() {}
        // set all values to uniform values.
        virtual void makeUniform() {}
                                                                                
        ///////////////////////////////////////////////////////////
        // read in the basic parameters, assuming file pointer
        // is located at the correct position.
        virtual void read(iDataStreamFile& is) {}
        ///////////////////////////////////////////////////////////
        // Do nothing.
        virtual void write(oDataStreamFile& os) {}
                                                                                
        void setLatticeADT(const LatticeADT &latticeAdt);
                                                                                
        ////////////////////////////////////////////////////////////////////////////
        // from base class EMable
        void emStartIteration() {}
        void emIncrement(logpr prob,vector < RV* >& parents, RV*r) {}
        void emEndIteration() {}
        void emSwapCurAndNew() {}
                                                                                
        // return the number of parameters for object.
        virtual unsigned totalNumberParameters() {return 0;}
                                                                                
        ///////////////////////////////////////////////////////////////
        // virtual functions for objects to do the actual work.
        virtual void emStoreObjectsAccumulators(oDataStreamFile& ofile) {}
        virtual void emLoadObjectsDummyAccumulators(iDataStreamFile& ifile) {}
        virtual void emZeroOutObjectsAccumulators() {}
        virtual void emLoadObjectsAccumulators(iDataStreamFile& ifile) {}
        virtual void emAccumulateObjectsAccumulators(iDataStreamFile& ifile) {}

protected:
        // returns the type of the sub-object in string
        // form that is suitable for printing and identifying
        // the type of the object.
        virtual const string typeName() {return std::string("LatticeEdgeCPT");}
                                                                                
        const LatticeADT *_latticeAdt;
};


#endif

