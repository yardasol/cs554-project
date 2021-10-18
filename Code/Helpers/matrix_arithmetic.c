#ifndef MATRIX_ARITHMETIC_C
#define MATRIX_ARITHMETIC_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "poisson_matrix_creaters.c" // to include A_csr and A_dia structs
// struct A_csr {
//     float* val;
//     float* col_ind;
//     float* row_ptr;
// };

// struct A_dia {
//     float** val;
//     float* off;
// };

/*Multipliers for three families of matrices*/
float* matVecProduct1(int n, float* y, float** x)
{
  float* b;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      b[j] += x[j][i] * y[i];
    }
  }
  return b;
}

float* matVecProductCSR1(int n, float* y, struct A_csr A)
{
  float* b; int col;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    b[i]=0;  
    for (int j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++)
    {
        col = A.col_ind[j];
        b[i] += A.val[j] * y[col];
    }
  }
  return b;
}


float* matVecProductDIA1(int n, int n_diag, float* y, struct A_dia A)
{
  float* b; 
  int offs, start, end;
  b = malloc(sizeof(float*) * n);
  for (int i=0; i<n; i++)
  {
    b[i]=0;  
    for (int j=0; j<n_diag; j++)
    {
        start = 0; end = n;
        offs = A.off[j];        
        if (offs>0){end -= offs;}
        else{start -= offs;}
        if((i>=start)&&(i<=end)){
            b[i] += A.val[j][i] * y[i+offs];
        }
    }
  }
  return b;
}

/*Inner product between two matrices*/
float innerProd1(int n, float* a, float* b)
{
  float c;
  for (int i=0; i<n; i++)
    c += a[i] * b[i];
  return c;      
}

/*Vector addition between matrices a and b, with a scale of k*/
float* VecAdd1(int n, float* a, float* b, float k)
{
  float* c;
  c = malloc(sizeof(float*) * n);    
  for (int i=0; i<n; i++)
    c[i] = a[i] + k*b[i];
  return c;      
}

#endif