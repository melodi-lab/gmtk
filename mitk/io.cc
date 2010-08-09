#include <iostream>
#include "mixNormalCollection.h"


///////////////////////  KMEANS READING PARAMETERS IN ///////////////

//////////////////// MixNormalCollection::readCurKMeansParamsBin ////////////////////

/**
 * read parameters from a file
 *
 * @param pf the input file pointer
 * @param forceAllActive whether set all GM to be active
 */
void MixNormalCollection::readCurKMeansParams(FILE *fp,bool isBin) {
  for ( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ) {
    if(isBin)
      p->readCurKMeansParamsBin(fp);
    else
      p->readCurKMeansParams(fp);
  }
}

//////////////////// MixNormal::readCurKMeansParamsBin ////////////////////

/**
 * read in the paramters to start in binary format
 *
 * @param fp the input binary file pointer
 * @exception NoSuchMethodException wrong file format
 * @throws ErrorReadingException error reading the file
 */
void MixNormal::readCurKMeansParamsBin(FILE *fp) {
  unsigned rc;

  unsigned dim = _numVariables;

  if ( (rc = fread(&_numVariables, sizeof(unsigned), 1, fp)) != 1 )
    error("Cannot read the number of variables from the input file");
  if ( (rc = fread(&_numMixtures, sizeof(unsigned), 1, fp)) != 1 )
    error("Cannot read the number of mixture components from input file");
  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    rc = fread(&(_alphas[l]), sizeof(PARAM_DATA_TYPE), 1, fp);
    rc += fread( (_means + l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    if ( rc != 1 + _numVariables )
      error("Cannot read alphas and means of mixture component # %d from input file",l);
    rc = 0;
    for ( unsigned i = 0; i < _numVariables; i++ )
      rc += fread( (_cov+l*dim*dim+i*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    if ( rc != _numVariables * _numVariables )
      error("Cannot read covs of mixture component # %d from input file",l);
      
  }
} // end readCurKMeansParamsBin


//////////////////// MixNormal::readCurKMeansParams ////////////////////

/**
 * read kmeans parameters an ascii file
 *
 * @param fp the input file pointer
 * @exception NoSuchMethodException wrong format of file
 * @throws ErrorReadingException error reading the file
 */
void MixNormal::readCurKMeansParams(FILE *fp) {
  unsigned tmp;
  int rc;
  unsigned dim = _numVariables;

  rc = fscanf(fp, "%u", &_numVariables);
  if ( rc == 0 || rc == EOF )
    error("cannot read the number of variables from parameter input file");

  rc = fscanf(fp, "%u", &_numMixtures);
  if ( rc == 0 || rc == EOF )
    error("cannot read the number of mixture components from parameter input file");

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    rc = fscanf(fp, "%u", &tmp);
    if ( rc == 0 || rc == EOF || tmp != l )
      error("cannot read the mixture component number %d from paramter input file",l);

    rc = fscanf(fp, "%le", &(_alphas[l]));
    if ( rc == 0 || rc == EOF )
      error("cannot read alpha value in mixture component # %d from parameter input file",l);

    for ( unsigned i = 0; i < _numVariables; i++ ) {
      rc = fscanf(fp, "%le",  (_means + l*dim+i));
      if ( rc == 0 || rc == EOF )
	error("cannot read means of mixture component # %d frm paramter input file",l);
    }
    
    for ( unsigned i = 0; i < _numVariables; i++ )
      for ( unsigned j = 0; j < _numVariables; j++ ) {
	rc = fscanf(fp, "%le",(_cov+l*dim*dim+i*dim+j));
	if ( rc == 0 || rc == EOF )
	  error("cannot read covs of mixture component # %d from paramter input file",l);
      }
  }
} // end readCurKMeansParams


//////  KMEANS WRITING PARAMETERS OUT ///////

void MixNormalCollection::writeCurKMeansParams(FILE *const fp,bool isBin){
  MixNormal *ftr_mi_p;    

  if (fseek (fp, 0L, SEEK_SET) != 0)
    error("Error seeking to beginning of output parameter file.");
  ftr_mi_p = _ftrMI;
  while (ftr_mi_p != _ftrMI_endp) {
    //      cout<<"Writing mixNormal KMeans parameters...\n";
      if(isBin)
	ftr_mi_p->printCurKMeansParamsBin(fp);
      else
	ftr_mi_p->printCurKMeansParams(fp);
      ftr_mi_p++;
  }
}


//////////////////// printCurKMeansParamsBin ////////////////////

/**
 * write to the c-like files in binary
 *
 * @param fp the output binary file pointer
 * @exception ErrorWritingException error writing to file
 */
void MixNormal::printCurKMeansParamsBin(FILE *fp) const {
  size_t rc;

  unsigned dim = _numVariables;

  if ( (rc = fwrite(&_numVariables, sizeof(unsigned), 1, fp)) != 1 )
    error("cannot write the number of variables to output file");
  if ( (rc = fwrite(&_numMixtures, sizeof(unsigned), 1, fp)) != 1 )
    error("cannot write the number of mixture components to output file");

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    rc = fwrite( (_alphas + l ), sizeof(PARAM_DATA_TYPE), 1, fp);
    rc += fwrite( (_means+l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    if ( rc != 1 + _numVariables )
      error("cannot write alphas and means to output file");

      for ( unsigned i = 0; i < _numVariables; i++ ) {
	rc = fwrite( (_cov+l*dim*dim+i*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
	if ( rc != _numVariables )
	  error("cannot write covs to output file");
      }

  }

} // end printCurKMeansParamsBin


//////////////////// printCurKMeansParams ////////////////////

/**
 * write kmeans parameters to an ascii file
 *
 * @param fp the output file pointer
 * @exception ErrorWritingException error writing to file
 */
void MixNormal::printCurKMeansParams(FILE *fp) const {
  unsigned dim = _numVariables;

  fprintf(fp, "%d\n", _numVariables);
  fprintf(fp, "%d\n", _numMixtures);

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    fprintf(fp, "%d\n", l);

    fprintf(fp, "%f\n", _alphas[l]);
    
    for ( unsigned i = 0; i < _numVariables; i++ )
      fprintf(fp, "%f ", *(_means+l*dim+i));
    fprintf(fp, "\n");

    for ( unsigned i = 0; i < _numVariables; i++ ) {
      for ( unsigned j = 0; j < _numVariables; j++ )
	fprintf(fp, "%f ", *(_cov+l*dim*dim+i*dim+j));
      fprintf(fp, "\n");
    }
  }
} // end printCurKMeansParams




///////////////// MixNormalCollection::WriteCurParams  /////////////////


/**
 * write the parameters EM out to a file in binary.  Only updated
 * (dirty) parameters are saved.  The others are skipped.  Information
 * about whether EM has converged for each given mixture is also
 * saved.
 * */

void MixNormalCollection::writeCurParams(FILE *const fp, bool isBin){
  if(fp == NULL) return;
  MixNormal *ftr_mi_p;    

  if (fseek (fp, 0L, SEEK_SET) != 0)
    error("Error seeking to beginning of output parameter file.");
  ftr_mi_p = _ftrMI;
  while (ftr_mi_p != _ftrMI_endp) {
    if(!isBin) 
      ftr_mi_p->printCurParams(fp);
    else {
      if (ftr_mi_p->dirty()) {
	ftr_mi_p->printCurParamsBin(fp);
	ftr_mi_p->reSetDirty();
      } else { // assume already saved.
	ftr_mi_p->seekOverCurParamsBin(fp);
      }
    }
    ftr_mi_p++;
  }

  if(!isBin) return;

  // save the active status as well
  ftr_mi_p = _ftrMI;
  while (ftr_mi_p != _ftrMI_endp) {
    const char act=1; const char inact=0;
    if (ftr_mi_p->active()) {
      if (fwrite(&act,sizeof(char),1,fp) != 1)
	error("Error writing active status.");
    } else {
      if (fwrite(&inact,sizeof(char),1,fp) != 1)
	error("Error writing active status.");
    }
    ftr_mi_p++;
  }
  if (fflush(fp) != 0)
    error("Error flushing mg parameters.");
}



//////////////////// MixNormal::printCurParamsBin ////////////////////

/**
 * write to the c-like files in binary
 *
 * @param fp the output binary file pointer
 * @exception ErrorWritingException error writing to file
 */
void MixNormal::printCurParamsBin(FILE *fp) const {
  unsigned i, l;
  size_t rc;

  unsigned dim = _numVariables;
  unsigned isFullCovar = (unsigned) _fullCoVar;

  //cout<<"_fullCoVar is "<<_fullCoVar<<" ("<<isFullCovar<<")\n"; 

  if ( (rc = fwrite(&_numVariables, sizeof(unsigned), 1, fp)) != 1 )
    error("cannot write to output file");
  if ( (rc = fwrite(&_numMixtures, sizeof(unsigned), 1, fp)) != 1 )
    error("cannot write to output file");
  if ( (rc = fwrite(&isFullCovar, sizeof(unsigned), 1, fp)) != 1 )
    error("cannot write to output file");

  for ( l = 0; l < _numMixtures; l++ ) {
    rc = fwrite( (_alphas + l ), sizeof(PARAM_DATA_TYPE), 1, fp);
    rc += fwrite( (_means+l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    rc += fwrite( (_invVars+l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);

    if ( rc != 1 + _numVariables + _numVariables )
      error("cannot write to output file");

    if ( _fullCoVar ) {
      for ( i = 0; i < _numVariables; i++ ) {
	rc = fwrite( (_b+l*dim*dim+i*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
	if ( rc != _numVariables )
	  error("cannot write to output file");
      }
    }
  }

  // seek over the remaining unused components
  if ( _numMixtures < _orgNumMixtures ) {
    unsigned offset = sizeof(PARAM_DATA_TYPE) * (1 + 2 * _numVariables);
    if ( _fullCoVar )
      offset += sizeof(PARAM_DATA_TYPE) * _numVariables * _numVariables;

    if ( fseek(fp, offset * (_orgNumMixtures - _numMixtures), SEEK_CUR) == -1 )
      error("problem seeking over mixture components");
  }
} // end writeCurParamsBin



//////////////////// printCurParams ////////////////////

/**
 * write EM parameters out in ASCII
 *
 * @param fp the output file pointer
 */
void MixNormal::printCurParams(FILE *fp) const {
  unsigned dim = _numVariables;

  fprintf(fp, "%d\n", _numVariables);
  fprintf(fp, "%d\n", _numMixtures);
  fprintf(fp, "%u\n", _fullCoVar);

  for ( unsigned l = 0; l < _numMixtures; l++ ) {
    fprintf(fp, "%d\n", l);
    fprintf(fp, "%f\n", _alphas[l]);

    for ( unsigned i = 0; i < _numVariables; i++ )
      fprintf(fp, "%f ", *(_means+l*dim+i));
    fprintf(fp, "\n");

    for ( unsigned i = 0; i < _numVariables; i++ )
      fprintf(fp, "%f ", *(_invVars+l*dim+i));
    fprintf(fp, "\n");

    if ( _fullCoVar ) {
      for ( unsigned i = 0; i < _numVariables; i++ ) {
	for ( unsigned j = 0; j < _numVariables; j++ )
	  fprintf(fp, "%f ", *(_b+l*dim*dim+i*dim+j));
	fprintf(fp, "\n");
      }
    }
  }
} // end MixNormal::printCurParams


//////////////////////////////// MixNormalCollection::ReadCurParams /////////////


/**
 * read parameters from a file
 *
 * @param pf the input file pointer
 * @param forceAllActive whether set all GM to be active
 */
int MixNormalCollection::readCurParams(FILE *fp, const bool forceAllActive, bool isBin) {
  //if(fp == NULL) return 0;
  int numActive = _numMis;

  for ( MixNormal *p = _ftrMI; p != _ftrMI_endp; p++ ) {
    if(isBin) 
      p->readCurParamsBin(fp);
    else
      p->readCurParams(fp);      
    p->reSetDirty();
  }

  if(!isBin) return numActive;

  // read the acitve status as well, if exists
  char active;
  if ( ! forceAllActive && (fread(&active, sizeof(char), 1, fp) == 1) ) {
    // then we presume there is status data
    numActive = 0;
    if ( active ) {
      _ftrMI->setActive();
      numActive++;
    } else
      _ftrMI->reSetActive();

    for ( MixNormal *p = _ftrMI + 1; p != _ftrMI_endp; p++ ) {
      if ( fread(&active, sizeof(char), 1, fp) != 1 )
	error("EOF encountered");
      if ( active ) {
	p->setActive();
	numActive++;
      } else
	p->reSetActive();
    }
  }

  return numActive;
} //end  MixNormalCollection::ReadCurParams


//////////////////// MixNormal::readCurParamsBin ////////////////////

/**
 * read in the paramters to start in binary format
 *
 * @param fp the input binary file pointer
 * @exception NoSuchMethodException wrong file format
 * @throws ErrorReadingException error reading the file
 */
void MixNormal::readCurParamsBin(FILE *fp) {
  unsigned i, l;
  unsigned rc;
  unsigned tmp;

  unsigned dim = _numVariables;
  DBGFPRINTF((stderr,"Before reading saved parameters;  _numVariables=%d and _numMixtures=%d (_orgNumMixtures=%d)\n",_numVariables,_numMixtures,_orgNumMixtures));
  if ( (rc = fread(&_numVariables, sizeof(unsigned), 1, fp)) != 1 )
    error("Cannot read the number of variables from the input file");
  if ( (rc = fread(&_numMixtures, sizeof(unsigned), 1, fp)) != 1 )
    error("Cannot read the number of mixture components from the input file");
  if ( (rc = fread(&tmp, sizeof(unsigned), 1, fp)) != 1 )
    error("Cannot read fullCoVar status from input file");
  DBGFPRINTF((stderr,"After reading saved parameters;  _numVariables=%d and _numMixtures=%d (_orgNumMixtures=%d)\n",_numVariables,_numMixtures,_orgNumMixtures));
  if ( tmp != (unsigned)_fullCoVar ) {
    warning("Conflicting settings of fullCoVar:  in mg file, fullCoVar is %d, but %d currently.",tmp,_fullCoVar);
  }
  for ( l = 0; l < _numMixtures; l++ ) {
    rc = fread(&(_alphas[l]), sizeof(PARAM_DATA_TYPE), 1, fp);
    rc += fread( (_means + l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    rc += fread( (_invVars + l*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
    if ( rc != 1 + _numVariables + _numVariables )
      error("Cannot read input file");

    if ( _fullCoVar ) {
      rc = 0;
      for ( i = 0; i < _numVariables; i++ )
	rc += fread( (_b+l*dim*dim+i*dim), sizeof(PARAM_DATA_TYPE), _numVariables, fp);
      if ( rc != _numVariables * _numVariables )
	error("Cannot read input file");
    }
  }

  // seek over the remaining unused components
  if ( _numMixtures < _orgNumMixtures ) {
    unsigned offset = sizeof(PARAM_DATA_TYPE) * (1 + 2 * _numVariables);
    if ( _fullCoVar )
      offset += sizeof(PARAM_DATA_TYPE) * _numVariables * _numVariables;

    if ( fseek(fp, offset * (_orgNumMixtures - _numMixtures), SEEK_CUR) == -1 )
      error("problem seeking over mixture components");
  }

  normalize();
} // end MixNormal::readCurParamsBin



//////////////////// MixNormal::readCurParams ////////////////////

/**
 * read the parameters in ascii
 *
 * @param fp the input file pointer
 * @exception NoSuchMethodException wrong format of file
 * @throws ErrorReadingException error reading the file
 */
void MixNormal::readCurParams(FILE *fp) {
  unsigned i, j, l;
  unsigned tmp;
  int rc;
  unsigned dim = _numVariables;

  rc = fscanf(fp, "%u", &_numVariables);
  if ( rc == 0 || rc == EOF )
    error("cannot read parameter input file");

  rc = fscanf(fp, "%u", &_numMixtures);
  if ( rc == 0 || rc == EOF )
    error("cannot read parameter input file");

  rc = fscanf(fp, "%u", &tmp);
  if ( rc == 0 || rc == EOF )
    error("cannot read parameter input file");
  if ( tmp != (unsigned)_fullCoVar )
    error("cannot change the property of gaussian mixture now");

  for ( l = 0; l < _numMixtures; l++ ) {
    rc = fscanf(fp, "%u", &tmp);
    if ( rc == 0 || rc == EOF || tmp != l )
      error("cannot read paramter input file");

    rc = fscanf(fp, "%le", &(_alphas[l]));
    if ( rc == 0 || rc == EOF )
      error("cannot read paramter input file");

    for ( i = 0; i < _numVariables; i++ ) {
      rc = fscanf(fp, "%le",  (_means + l*dim+i));
      if ( rc == 0 || rc == EOF )
	error("cannot read paramter input file");
    }

    for ( i = 0; i < _numVariables; i++ ) {
      rc = fscanf(fp, "%le", (_invVars + l*dim+i));
      if ( rc == 0 || rc == EOF )
	error("cannot read paramter input file");
    }

    if ( _fullCoVar ) {
      for ( i = 0; i < _numVariables; i++ )
	for ( j = 0; j < _numVariables; j++ ) {
	  rc = fscanf(fp, "%le",(_b+l*dim*dim+i*dim+j));
	  if ( rc == 0 || rc == EOF )
	    error("cannot read paramter input file");
	}
    }
  }

  normalize();
} // end MixNormal::readCurParams



//////////////////// MixNormal::seekOverCurParamsBin ////////////////////

/**
 * seek over the file of parameters
 *
 * @param fp the file pointer to parameters
 * @throws IOException wrong size of the file
 */
void MixNormal::seekOverCurParamsBin(FILE *fp) const {
  unsigned offset = sizeof(PARAM_DATA_TYPE) * (1 + 2 * _numVariables);
  if ( _fullCoVar )
    offset += sizeof(PARAM_DATA_TYPE) * _numVariables * _numVariables;

  if ( fseek(fp, sizeof(unsigned) * 3 + offset * _orgNumMixtures,
	     SEEK_CUR) == -1 )
    error("problem seeking over mixture components in file");
} // end MixNormal::seekOverCurParamsBin




//////////////////// print ////////////////////

/**
 *  print current EM parameters
 * 
 */
void MixNormal::print() {

  unsigned l, i, j;
  unsigned dim = _numVariables;
  for( l = 0; l < _numMixtures; l++){
    cout << "Mixture " << l << endl;
    cout << "weight " << _alphas[l] << endl;
    cout << "mean ";
    for( i = 0; i < _numVariables; i++) cout << *(_means+l*dim+i) << " ";
    cout << endl;
    cout << "invVars ";
    for( i = 0; i < _numVariables; i++) cout << *(_invVars+l*dim+i) << " ";
    cout << endl;

    if( _fullCoVar ) {
      cout << "b ";
      for( i = 0; i < _numVariables; i++){
	for( j = 0; j < _numVariables; j++) cout << *(_b + l*dim*dim + i*dim + j) << " ";
	cout << endl;
      }
    }
    cout << endl;
  }
}


//////////////////// printParams ////////////////////

/**
 * print the paramters for debugging only
 */
void MixNormal::printParams() const {
  cout << "alphas" << endl;
  for ( unsigned l = 0; l < _numMixtures; l++ )
    cout << _alphas[l] << " ";
  cout << endl;

  for ( unsigned n = 0; n < _numVariables; n++ ) {
    cout << "variable" << n << endl;		// output the variable index

    cout << "mixture\t mean\t variance" << endl;
    for ( unsigned l = 0; l < _numMixtures; l++ ) {
      cout << l << "\t" <<  *(_means + l*_numVariables+n) << "\t" << *(_invVars +l*_numVariables+n);
    }
    cout << endl;
  }
} // end printParams


/**
 * when a matrix is not positive definite, dumps relevant data to a file
 *
 */

void MixNormal::notPosDef(unsigned index, unsigned compNum, unsigned dim, PARAM_DATA_TYPE* cov) {
  FILE* ofp = fopen(ERR_FILE,"w");
  if(ofp != NULL) {
    fprintf(ofp,"ERROR: in MixNormal::endEpoch(),  tuple # %d, mixture component # %d, the covariance matrix is singular.\nCov [%d x %d] =", index, compNum,dim,dim);
    for ( unsigned i = 0; i < dim*dim; i++ ) 
      fprintf(ofp,"%f ",*(cov+i));
    fprintf(ofp,"\n");
  }
  error("In MixNormal::endEpoch(), the covariance matrix is singular.  Error log written to %s",ERR_FILE);
} 
