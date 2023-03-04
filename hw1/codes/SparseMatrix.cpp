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
            sparse_matrix[i][j]=0;
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

// accessor 
int SparseMatrix::GetSize(int i) const
{
    assert(i==1 || i==2); // only two values of i are accepted 
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

// overload assignment operator
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& otherSparseMatrix)
{
    this->m = otherSparseMatrix.m;
    this->n = otherSparseMatrix.n; 
    this->nnz = otherSparseMatrix.nnz; 
    this->colptr = otherSparseMatrix.colptr;
    this->rowidx = otherSparseMatrix.rowidx;
    this->nzval = otherSparseMatrix.nzval;

    return *this;
}

// unary +
SparseMatrix SparseMatrix::operator+() const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
    {
        new_nzvalues[i] = + this->nzval[i]; 
    }

    SparseMatrix new_plus_sm(this->nnz, this->rowidx, this->colptr, new_nzvalues, this->GetSize(1), this->GetSize(2));
    return new_plus_sm;
}

// unary -
SparseMatrix SparseMatrix::operator-() const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
    {
        new_nzvalues[i] = - this->nzval[i]; 
    }

    SparseMatrix new_minus_sm(this->nnz, this->rowidx, this->colptr, new_nzvalues, this->GetSize(1), this->GetSize(2));
    return new_minus_sm;
}

// scalar multiplication
SparseMatrix SparseMatrix::operator *(double a) const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
    {
        new_nzvalues[i] = a * this->nzval[i]; 
    }

    SparseMatrix new_multiply_sm(this->nnz, this->rowidx, this->colptr, new_nzvalues, this->GetSize(1), this->GetSize(2));
    return new_multiply_sm;
}