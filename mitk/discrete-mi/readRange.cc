#include "readRange.h"

#define MAX_QUEUE_SIZE 100

int parseSubExpr(char* &lp, int queue[], int &qInd);


//creates the two set of features and their corresponding lags, of which we want to compute the MI 
int RangeSet::createSets(char* line) {
  
  int status;
  char *lp=line;
  char *pos;

  if( ( pos = strchr(lp,'#') ) != NULL ) *pos = '\0';  //remove comments
  while(isspace(*lp)) ++lp;  //eat spaces
  if(*lp=='\0' || *lp =='\n') return COMMENT;
  if(*lp != '[') {
    //printf("ERROR: unrecongnized start token ('%c') for set A\n",*lp); 
    return UNRECOGNIZED_TOKEN; 
  }

  //check [ ] [ ] format
  pos=lp;
  for(int i=0;i<2;++i) {
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

  if(*lp != '[') {
    //printf("ERROR: unrecongnized start token  ('%c' (%d ascii)) for set B ( rest of line = %s)\n",*lp,*lp,lp); 
    return UNRECOGNIZED_TOKEN;
  }

  if( (status = createOneSet(1,lp) ) != SUCCESS) return status; //B = set 1

  if(size[0]==0 || size[1]==0) return EMPTY_SET; 
  
  //adjust the lag so that the smallest is >=0
  for(int i=0;i<size[1];++i) {
    if ( set[1][i].lag > max_lag )
      max_lag = set[1][i].lag;
    if ( set[1][i].lag < min_lag )
      min_lag = set[1][i].lag;
  }

  //shift everything so that the smallest lag is 0
  if(min_lag < 0) {
    for( int i=0; i < size[0]; ++i)  set[0][i].lag += (-min_lag);  //A
    for( int i=0; i < size[1]; ++i)  set[1][i].lag += (-min_lag);  //B
    max_lag += (-min_lag);
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
    if(fqInd <1 ||  lqInd <1) return EMPTY_SET;
    if(lqInd == 1) {
      lqInd = fqInd;
      for(int i=1;i<lqInd;++i) lagQueue[i] = lagQueue[0];
    }
    else if(lqInd != fqInd) {
      printf("Number of features does not match number of lags.\n");
      exit(-1);
    }
    //check redundancy
    for(int i=0; i< getSize(setNum);++i) {
      for(int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	if(set[setNum][i].feat == featQueue[fqInd2] && set[setNum][i].lag == lagQueue[lqInd2]) {
	  printf("ERROR: redundant feature/lag pair\n");
	  exit(-1);
	}
      }
    }

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
    //store into set[setNum].  Could also do it while checking for redundancy
    for( int fqInd2=0,lqInd2=0; fqInd2 < fqInd; ++fqInd2, ++lqInd2) {
	set[setNum][size[setNum]].feat = featQueue[fqInd2];
	set[setNum][size[setNum]].lag = lagQueue[lqInd2];
	++size[setNum];
    }


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
    //printf("invalid range format: %s\n", lp);
    return INVALID_FORMAT;
  }
  return SUCCESS;
}
