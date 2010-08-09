#include "readRange.h"

#define MAX_QUEUE_SIZE 100

int parseSubExpr(char* &lp, int queue[], int &qInd);


//creates the two set of features and their corresponding lags, of
//which we want to compute the MI
int RangeSet::createSets(char* line) {
  
  int status;
  char *lp=line;
  char *pos;

  if( ( pos = strchr(lp,'#') ) != NULL ) *pos = '\0';  //remove comments
  while(isspace(*lp)) ++lp;  //eat spaces
  if(*lp=='\0' || *lp =='\n') return COMMENT;
  if(*lp != '[') {
    return UNRECOGNIZED_TOKEN; 
  }

  //check [ ] [ ] format
  pos=lp;
  for(int i=0;i<1;++i) {  // changed i<2 to i<1.  We want to allow [ ] tuples only
    if( ( pos = strchr(pos,'[') ) == NULL )
      return BRACKET_ERROR;
    if( ( pos = strchr(pos,']') ) == NULL )
      return BRACKET_ERROR;
  }

  if( (status = createOneSet(0,lp) ) != SUCCESS) return status;  //A = set 0
  


  for(int i=0;i<size[0];++i) {
    if ( set[0][i].lag > max_lag )
      max_lag = set[0][i].lag;
    if ( set[0][i].lag < min_lag )
      min_lag = set[0][i].lag;
  }

  while( isspace(*lp) ) ++lp;  //eat spaces

  //#if ALLOW_ENTROPY
  if( ( pos = strchr(pos,'[') ) == NULL ) {
    // we have a [ ] only format
    size[1]=0;
    return SUCCESS;
  }
  //#endif


  if(*lp != '[') {
    return UNRECOGNIZED_TOKEN;
  }

  if( (status = createOneSet(1,lp) ) != SUCCESS) return status; //B = set 1

  if(size[0]==0 || size[1]==0) { 
    if(size[0]==0)
      DBGFPRINTF((stderr,"set 0 is empty.\n"));
    if(size[1]==0)
      DBGFPRINTF((stderr,"set 1 is empty.\n"));
    return EMPTY_SET; 
  }
  
  for(int i=0;i<size[1];++i) {
    if ( set[1][i].lag > max_lag )
      max_lag = set[1][i].lag;
    if ( set[1][i].lag < min_lag )
      min_lag = set[1][i].lag;
  }

  return SUCCESS;
}


//create one of the two sets above
int RangeSet::createOneSet(int setNum, char* &lp) {

  int featQueue[MAX_QUEUE_SIZE], lagQueue[MAX_QUEUE_SIZE];
  int fqInd = 0, lqInd = 0;
  int status;

  while(*lp != ']') {
    ++lp;
    
    fqInd=lqInd=0;
    while(*lp != '@') {
      while( isspace(*lp) ) ++lp;
      //queue the features until you read the lags
      status = parseSubExpr(lp, featQueue,fqInd);
      if(status != SUCCESS) return status;
      ++lp;
    }

    ++lp;
    while(*lp !=';' && *lp != ']') {
      while( isspace(*lp) ) ++lp;  //eat spaces
      //queue the lags
      status = parseSubExpr(lp, lagQueue,lqInd);
      if(status != SUCCESS) return status;
      ++lp;
    }
    
    //check #features == #lags (special case #lags == 1 => repeat it for all features)
    if(fqInd <1 ||  lqInd <1) {
      DBGFPRINTF((stderr,"Found empty set.\n"));
      return EMPTY_SET;
    }
    if(lqInd == 1) {
      lqInd = fqInd;
      for(int i=1;i<lqInd;++i) lagQueue[i] = lagQueue[0];
    }
    else if(lqInd != fqInd) {
      printf("Number of features does not match number of lags.\n");
      exit(-1);
    }
    //check redundancy
    int set_size_so_far=getSize(setNum);
    if(set_size_so_far==0) {
      for( int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	set[setNum][size[setNum]].feat = featQueue[fqInd2];
	set[setNum][size[setNum]].lag  = lagQueue[lqInd2];
	++size[setNum];
      }
    }
    else {
      for(int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	int i;
	for(i=0; i< set_size_so_far;++i) {
	  if(set[setNum][i].feat == featQueue[fqInd2] && set[setNum][i].lag == lagQueue[lqInd2]) {
#if ALLOW_REDUNDANT_PAIRS
	    printf("WARNING: redundant feature/lag pair (%d@%d).  Will ignore it.\n",set[setNum][i].feat, set[setNum][i].lag);
	    break;
#else
	    printf("ERROR: redundant feature/lag pair\n");
	    exit(-1);
#endif
	  }
	}
	if(i==set_size_so_far) { // we did not find duplicates
	    set[setNum][size[setNum]].feat = featQueue[fqInd2];
	    set[setNum][size[setNum]].lag  = lagQueue[lqInd2];
	    ++size[setNum];
	  }
	
      }
    }
#if CHECK_TUPLE_OVERLAP
    //check overlap
    if(setNum==1) { //no need to do it for the first set since the second one is empty then
      for(int i=0; i< getSize((setNum+1)%2);++i) {
	for( int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	  if(set[(setNum+1)%2][i].feat == featQueue[fqInd2] && set[(setNum+1)%2][i].lag == lagQueue[lqInd2]) {
	    printf("ERROR: feature/lag pair overlap\n");
	    exit(-1);
	  }
	}
      }
    }
#endif
    //store into set[setNum].  Could also do it while checking for redundancy
    // It's now done while checking for redundancy
    /*
    for( int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	set[setNum][size[setNum]].feat = featQueue[fqInd2];
	set[setNum][size[setNum]].lag  = lagQueue[lqInd2];
	++size[setNum];
    }
    */

  }
  if(fqInd==0) return EMPTY_SET;
  lp++;
  return SUCCESS;
}

int parseSubExpr(char* &lp, int queue[], int &qInd) {
  int scanRead, num,step;

  if(isdigit(*lp) || (*lp == '-' && isdigit(*(lp+1))) ) {
    //read first number m 
    scanRead = sscanf(lp,"%d",&num);
    assert(scanRead==1);
    while(isdigit(*lp) || (*lp == '-' && isdigit(*(lp+1)))) ++lp;
    while( isspace(*lp) ) ++lp;
    queue[qInd++] = num;
    assert(qInd < MAX_QUEUE_SIZE);
    
    if(*lp==':') {  //range m:s:n or m:n
      lp++;
      while( isspace(*lp) ) ++lp;
      if( !isdigit(*lp) && !(*lp == '-' && isdigit(*(lp+1))) ) return UNRECOGNIZED_TOKEN;
      
      //get step, s, or end of range, n, if no step is specified
      scanRead = sscanf(lp,"%d",&num);
      assert(scanRead==1);
      while(isdigit(*lp) || (*lp == '-' && isdigit(*(lp+1)))) ++lp;
      while( isspace(*lp) ) ++lp;
      
      //got one num.  Get second one if any
      if(*lp==':') { //last number was the step.  Get the end of range.
	lp++;
	while( isspace(*lp) ) ++lp;
	if( !isdigit(*lp) && !(*lp == '-' && isdigit(*(lp+1))) ) return UNRECOGNIZED_TOKEN;
	step = num;
	scanRead = sscanf(lp,"%d",&num);
	assert(scanRead==1);
	while(isdigit(*lp) || (*lp == '-' && isdigit(*(lp+1)))) ++lp;
	while( isspace(*lp) ) ++lp;
      }
      else  //last number was end of range
	step=1;
          
      //queue  m+1 through n
      int lastInd = qInd-1;
      for(int i=queue[lastInd]+step; i<= num;i+=step) {
	queue[qInd++] = i;
	assert(qInd < MAX_QUEUE_SIZE);
      }
    }
    --lp;  //will be incremented in the caller
  }
  else if(*lp !=',') {
    return INVALID_FORMAT;
  }
  return SUCCESS;
}
