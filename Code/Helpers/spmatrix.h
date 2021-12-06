#ifndef SPMATRIX_H
#define SPMATRIX_H

struct A_csr {
    float* val;
    int* col_ind;
    int* row_ptr;
};

struct A_coo_entry {
    float val;
    unsigned int row;
    unsigned int col;
};

struct A_coo {
    struct A_coo_entry* data;
    unsigned int nnz;
    unsigned int rows;
    unsigned int cols;
};

struct A_dia {
    float** val;
    float* off;
};

#endif /* SPMATRIX_H */
