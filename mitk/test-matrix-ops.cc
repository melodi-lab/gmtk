#include "matrix-ops.h"
//#include "data-points.h"

bool testLDU(double* cov, unsigned dim) {
  //moment to LDU
  double *newCov =new double[dim*dim];
  double *b =new double[dim*dim];
  double *tmp =new double[dim*dim];
  double *invVars = new double [dim];

  cout<<"BEFORE DIAG_CHOL\n";
  cout<<"InvVars: ";
  for ( unsigned i = 0; i < dim; i++) cout<<*(invVars+i)<<"  ";
  cout<<endl;
  cout<<"B: ";
  for ( unsigned i = 0; i < dim*dim; i++) cout<<*(b+i)<<"  ";
  cout<<endl;
  cout<<"cov: ";
  for ( unsigned i = 0; i < dim*dim; i++) cout<<*(cov+i)<<"  ";
  cout<<endl;

  diagCholeskyDecomp(cov,b,invVars,dim);

  

  //  exit(0);

  invertVec(invVars,invVars,dim);
  inverse(b,tmp,dim); 
  transpose(b,tmp,dim);


  //inverse(b,b,dim); 
  //transpose(tmp,b,dim);

  
  // LDU to moment
  //inverse(tmp,b,dim);
  inverse(b,b,dim);
  for ( unsigned i = 0; i < dim; i++) 
    for ( unsigned j = 0; j < dim; j++) 
      for ( unsigned k = 0; k< dim; k++) 
	*(newCov + i*dim + j) += *(b + i*dim + k) / *(invVars+k) *  *(b + j*dim + k);

  for ( unsigned i = 0; i < dim*dim; i++) 
    cout << *(cov+i) <<"\t"<<*(newCov+i)<<"\t"<<*(b+i)<<"\t"<<*(tmp+i)<<endl;
  
  return true;
}

bool testTranspose(double * M,unsigned dim) {

  double* N = new double[dim*dim];

  transpose(N,M,dim);
  transpose(N,N,dim); 
  
  return true;

}

int main(void) {
  double * M = new double [4];
  double * N = new double [4];


  M[0]= 6;
  M[1] = 2;
  M[2] = 2;
  M[3] = 4;

  N[0]=1;
  N[1]=2;
  N[2]=3;
  N[3]=4;

  testLDU(M,2);

  return 0;
}
