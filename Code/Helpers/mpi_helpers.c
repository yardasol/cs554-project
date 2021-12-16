#ifndef MPI_HELPERS_C
#define MPI_HELPERS_C

#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "poisson_matrix_creaters.c"
#include "helpers.c"

struct par_multdat {
    int* offsetv_d ; /* used to determine start row sent to each worker */
    int rows_d ; /* no. of rows of matrix A sent to each worker */
    int rank_d ;
    int n_d ;
    int nwrks_d ;
    int n_diag_d;
};

struct par_multdat parmult_struct_assign(int *offsvec, int rows, int rank, int N, int nwrks, int n_diag)
{
    struct par_multdat pmdat;
    pmdat.offsetv_d = offsvec;
    pmdat.rows_d = rows;
    pmdat.rank_d = rank;
    pmdat.n_d = N;
    pmdat.nwrks_d = nwrks;
    pmdat.n_diag_d = n_diag;
    return pmdat;
}

struct par_multdat parmult_create(int N) {
    int numprocs, rank, rc;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    int nwrks, rows;
    nwrks = numprocs-1; // no. of workers

    int* offsvec =  row_load_allot(N,numprocs); // vector of offsets for each worker
    if (rank>0) { rows = offsvec[rank]-offsvec[rank-1]; } //Calc. no. of rows of A this worker rank deals with

    return parmult_struct_assign(offsvec, rows, rank, N, nwrks, 0);
}


/*MPI Multipliers for 3 matrix families*/
float* mpiMatVecProduct1( struct par_multdat pmd, float* x, float** a ) {

    int	source, dest, rows, offset, i, j, ind, ioff;

    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    float* b = create1dZeroVec(n);
    MPI_Status status;

    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = offsv[i-1]; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&b[ind], rows, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(b, n, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return b;
    }

    if (rank > 0) {
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            for (j=0; j< n; j++)
            { b[ioff] = b[ioff] + a[ioff][j] * x[j]; }
        }
        MPI_Send(&b[offset], rows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(b, n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return b;
    }
}


float* mpiMatVecProductCSR1( struct par_multdat pmd, float* y, struct A_csr A)
{
    int source, dest, rows, offset, i, j, ind, ioff, col;

    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    float* b = create1dZeroVec(n);
    MPI_Status status;

    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = offsv[i-1]; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&b[ind], rows, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(b, n, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return b;
    }

    if (rank > 0) {
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            for (int j=A.row_ptr[ioff]; j<A.row_ptr[ioff+1]; j++)
            {
                col = A.col_ind[j];
                b[ioff] += A.val[j] * y[col];
            }
        }
        MPI_Send(&b[offset], rows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(b, n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return b;
    }
}

float* mpiMatVecInvDProductCSR1( struct par_multdat pmd, float* y, struct A_csr A)
{
    int source, dest, rows, offset, i, j, ind, ioff, col;

    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    float* b = create1dZeroVec(n);
    MPI_Status status;
    float diag;

    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = offsv[i-1]; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&b[ind], rows, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(b, n, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return b;
    }

    if (rank > 0) {
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            diag = 0.f;
            for (int j=A.row_ptr[ioff]; j<A.row_ptr[ioff+1]; j++) {
                if (A.col_ind[j] == ioff) {
                    diag = A.val[j];
                    break;
                }
            }

            for (int j=A.row_ptr[ioff]; j<A.row_ptr[ioff+1]; j++)
            {
                col = A.col_ind[j];
                b[ioff] += (A.val[j] / diag) * y[col];
            }
        }
        MPI_Send(&b[offset], rows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(b, n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return b;
    }
}

float* mpiMatVecProductDIA1( struct par_multdat pmd, float* y, struct A_dia A)
{
    int   source, dest, rows, offset, i, j, ind, ioff;
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    int n_diag = pmd.n_diag_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    float* b = create1dZeroVec(n);
    MPI_Status status;

    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = offsv[i-1]; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&b[ind], rows, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(b, n, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return b;
    }

    if (rank > 0) {
        int offs, start, end;
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            for (j=0; j< n_diag; j++){
                start = 0; end = n;
                offs = A.off[j];
                if (offs>0){end -= offs;}
                else{start -= offs;}
                if((ioff>=start)&&(ioff<=end)){
                    b[ioff] +=  A.val[j][ioff] * y[ioff+offs];
                }
            }
        }
        MPI_Send(&b[offset], rows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(b, n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return b;
    }
}

/*MISC mpi arithmetic*/
float* mpiVecAdd1( struct par_multdat pmd, float* a, float* b, float k)
{
    int source, dest, rows, offset, i, j, ind, ioff;
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    MPI_Status status;

    float* c = create1dZeroVec(n);

    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = offsv[i-1]; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&c[ind], rows, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(c, n, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return c;
    }

    if (rank > 0) {
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            c[ioff] = a[ioff] + k*b[ioff];
        }
        MPI_Send(&c[offset], rows, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(c, n, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return c;
    }
}

float mpiInnerProd1( struct par_multdat pmd, float* a, float* b)
{
    int source, dest, rows, offset, i, j, ind, ioff;
    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    MPI_Status status;

    float c = 0;
    float* carr = create1dZeroVec(numworkers);
    if (rank == 0)
    {
        for (i=1; i<=numworkers; i++) {
            source = i; ind = i-1; rows = offsv[i]-offsv[i-1];
            MPI_Recv(&carr[ind], 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
            c += carr[ind];
        }
        for (i=1; i<=numworkers; i++) {MPI_Send(&c, 1, MPI_FLOAT, i, 2, MPI_COMM_WORLD);}
        return c;
    }

    if (rank > 0) {
        for (i = 0; i<rows; i++)
        {
            ioff = i+offset;
            c += a[ioff] * b[ioff];
        }
        MPI_Send(&c, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(&c, 1, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
        return c;
    }
}


/*MPI Algorithms*/
float* mpiTriangularSolveCSR1( struct par_multdat pmd, float* r, struct A_csr L, struct A_csr U)
{
    int source, dest, rows, offset, i, j, ind, ioff, col;

    int rank = pmd.rank_d;
    int n = pmd.n_d;
    int numworkers = pmd.nwrks_d;
    int* offsv = pmd.offsetv_d;
    offset = pmd.offsetv_d[rank-1];
    rows = pmd.rows_d;
    float* y = create1dZeroVec(n);
    float* z = create1dZeroVec(n);
    MPI_Status status;

    if (rank == 0)
    {
        int row, col;
        int row_ptr, next_row_ptr;
        int i;
        // Ly = r
        for (row = 0; row < n; row++){
            row_ptr = L.row_ptr[row];
            next_row_ptr = L.row_ptr[row+1];
            for (i = row_ptr; i < next_row_ptr-1; i++){
                col = L.col_ind[i];
                y[row] -= L.val[i] * y[col];
            }
            y[row] += r[row];
            y[row] = y[row] / L.val[i];
        }

        //Uz = y
        for (row = n-1; row >= 0; row--){
            row_ptr = U.row_ptr[row+1];
            next_row_ptr = U.row_ptr[row];
            for (i = row_ptr-1; i > next_row_ptr; i--){
                col = U.col_ind[i];

                z[row] -= U.val[i] * z[col];
            }
            z[row] += y[row];
            z[row] = z[row] / U.val[i];
        }

    }

    MPI_Bcast(z, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    return z;
}

/*MPI Verifiers*/
void mult_Output_verify(int N, int n_diag, int numprocs, int rank, int* offsvec,
        float** A, struct A_csr A_CSR, struct A_dia A_DIA, float* x, float*b,
        float* b_csr, float* b_dia, float* b_p, float* b_p_csr, float* b_p_dia)
{
    if (rank == 0)
    {
        //Output
        print_mat(N,N,A);
        print_vec(N,x);
        print_vec(N,b_p);
        print_vec(N,b_p_csr);
        /*print_vec(N,b_p_dia);
          print_vec(N,b); print_vec(N,b0);
          print_vec(N,b_csr); print_vec(N,b_dia);*/
        b = matVecProduct1(N,x,A);
        b_csr = matVecProductCSR1(N,x,A_CSR);
        b_dia = matVecProductDIA1(N,n_diag,x,A_DIA);
        sermult_debugger(N,b_csr,b); sermult_debugger(N,b_dia,b);
        parmult_debugger(N,numprocs,b_p,b,offsvec);
        parmult_debugger(N,numprocs,b_p_csr,b,offsvec);
        parmult_debugger(N,numprocs,b_p_dia,b,offsvec);
    }
    return;
}



#endif
