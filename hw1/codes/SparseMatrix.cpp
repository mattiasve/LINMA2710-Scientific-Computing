#include "SparseMatrix.hpp"
#include <cassert>
#include <iostream>

// zero sparse matrix contructor
SparseMatrix::SparseMatrix(int numRows, int numCols)
{
    int sparse_matrix [numRows][numCols]; 
    for (int i=0; i<numRows; i++)
    {
        for (int j=0; j<numCols; j++)
            {sparse_matrix[i][j]=0; }
    }
}

// deep copy constructor
SparseMatrix::SparseMatrix(int m, int n, int nnz, int* colptr, int* rowidx, double* nzval)
{
    this->m = m;
    this->n = n;
    this->nnz = nnz; 
    this->colptr = new int[n+1]; 
    this->rowidx = new int[nnz]; 
    this->nzval = new double [nnz]; 

    // copy date inside arrays
    for (int i=0; i<n+1; i++)
        this->colptr[i] = colptr[i]; 
    
    for (int j=0; j<nnz; j++)
    {
        this->rowidx[j] = rowidx[j]; 
        this->nzval[j] = nzval[j]; 
    }
}

// specialized contructor
SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int size1, int size2) 
{
    return 0.0;
}

// accesor 
int SparseMatrix::GetSize(int i) const
{
    assert(i==0 || i==1); // only two values of i are accepted 
    if (i == 1)
        return this->m;
    else
        return this->n;
}

// destructor
SparseMatrix::~SparseMatrix()
{
    delete[] colptr; 
    delete[] rowidx; 
    delete[] nzval; 
}