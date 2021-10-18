#ifndef HELPERS_C
#define HELPERS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


float VecErrNorm(int n, float* a, float* b)
{
  // Compare the L-infinity norm of sol a w.r.t. exact b
  float tmp, errnorm; 
  float max = 0; 
  float maxb = 0;
  for (int i=0; i<n; i++)
  {
    tmp = fabs(a[i] - b[i]); 
    max = (tmp>max) ? tmp : max;
    maxb = (fabs(b[i])>maxb) ? fabs(b[i]) : maxb;
  }
  errnorm = max/maxb;
  //printf("%3f  , %3f , %3f", max,maxb,errnorm);printf("\n  ");
  return errnorm;      
}

void print_mat(int nr, int nc, float** m)
{
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<nc; j++)
      printf("%3f  ", m[i][j]);
    printf("\n");
  }
  printf("\n");
}

void print_vec(int n, float* v)
{
  for (int i=0; i<n; i++)
    printf("%3f \n", v[i]);
  printf("\n");
}

int* row_load_allot(int n, int ptot)
{
    int nwrks, offset, avrow, rows, lrow, roweq, rowrem;
    int* offsv = malloc(sizeof(int*) * ptot);    

    nwrks = ptot-1; // no. of workers
    avrow = floor((float)n/(float)nwrks); // average no. of rows of A each worker deals with      
    lrow = 0;
    roweq = avrow*nwrks;
    if (roweq<n){lrow=1; rowrem = n-roweq;} 
    offset = 0; offsv[0]= 0;

    for (int k=1; k<ptot; k++)
    {
        if(k>rowrem){lrow=0;}
        rows = avrow + lrow; 
        offset = offset + rows; 
        offsv[k] = offset; 
    }

    return offsv;
}


void communicate_xvec(int N, int rank, int nwrks, float* x)
{
    int source, dest, tag;
    MPI_Status status; MPI_Request request;
    if (rank == 0)
    {  
        tag = 1;
        MPI_Request req[nwrks]; MPI_Status stat[nwrks];
        for (dest=1; dest<=nwrks; dest++)
            { MPI_Isend(&x[0], N, MPI_FLOAT, dest, tag, MPI_COMM_WORLD, &req[dest-1]);}

        MPI_Waitall(nwrks, req, stat);
    } else if (rank>0)
    {
        source = 0; tag = 1;
        MPI_Recv(&x[0], N, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    }
    //if (rank>0){MPI_Wait(&request, MPI_STATUS_IGNORE);}

    MPI_Barrier(MPI_COMM_WORLD); // make sure everybody has received the same x by here
    return;
}

void sermult_debugger( int n, float* b, float* b0 ) {
    float multerr = VecErrNorm(n,b,b0);
    printf("Error b/w the diff serial mults: %f \n \n", multerr);
    return;
}

void parmult_debugger(int n, int numprocs, float* b, float* b0, int* offsv) {
    float multerr = VecErrNorm(n,b,b0);
    for(int ii=0;ii<numprocs;ii++){printf("offsvec[%d]=%d\n",ii,offsv[ii]);}
    printf("Error b/w serial and parallel mults: %f \n \n", multerr);
    return;
}


#endif