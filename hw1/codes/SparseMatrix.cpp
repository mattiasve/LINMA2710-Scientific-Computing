#include "SparseMatrix.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

// zero sparse matrix contructor
SparseMatrix::SparseMatrix(int numRows, int numCols)
{
    m = numRows; 
    n = numCols; 
    nnz = 0; 
    colptr = NULL;
    rowidx = NULL;
    nzval = NULL;
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

    // copy data inside arrays
    for (int i=0; i<n+1; i++)
        this->colptr[i] = colptr[i]; 
    
    for (int j=0; j<nnz; j++)
    {
        this->rowidx[j] = rowidx[j]; 
        this->nzval[j] = nzval[j]; 
    }
}

// function used to sort arrays in the constructor
void selectionSort(int* to_sort_by, int* sorted, double* nvals, int low, int high) 
{
    int i, j, idx;

    for (i = low; i < high - 1; i++) {
        idx = i;
        for (j = i + 1; j < high; j++)
            if (to_sort_by[j] < to_sort_by[idx])
                idx = j;

        std::swap(to_sort_by[idx], to_sort_by[i]);
        std::swap(sorted[idx], sorted[i]);
        std::swap(nvals[idx], nvals[i]);
    }
}

// constructor
SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzvalues, int size1, int size2)
{
    this->nnz = nnz;
    this->m = size1 ; 
    n = size2 ;
    colptr = new int[n+1]; //leak
    rowidx = new int[nnz]; //leak
    nzval  = new double[nnz]; //leak
    int colidx[nnz];

    std::copy(cidx, cidx + nnz, colidx);
    std::copy(ridx, ridx + nnz, rowidx);
    std::copy(nzvalues, nzvalues + nnz, nzval);

    selectionSort(colidx, rowidx, nzval, 0, nnz);
    int j = 0;
    int i = 0;
    int current_col = 0;
    *(colptr + current_col) = i;
    while (i<nnz-1)
    {
        j = i;
        *(colptr+current_col) = i;
        while (*(colidx + i) == *(colidx+ j + 1)) {j++;}
        selectionSort(rowidx, colidx, nzval, i, j);
        i = j + 1;
        current_col++;
    }
    *(colptr + current_col) = i;
    *(colptr + n) = nnz;
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
    assert(this->m == otherSparseMatrix.GetSize(1)
    && this->n == otherSparseMatrix.GetSize(2));

    this->m = otherSparseMatrix.m;
    this->n = otherSparseMatrix.n; 
    this->nnz = otherSparseMatrix.nnz; 
    this->colptr = new int[otherSparseMatrix.n+1];
    this->rowidx = new int[otherSparseMatrix.nnz];
    this->nzval = new double[otherSparseMatrix.nnz];

    for (int i=0; i<otherSparseMatrix.nnz; i++)
    {
        *(rowidx+i) = otherSparseMatrix.rowidx[i];
        *(nzval+i) = otherSparseMatrix.nzval[i];
    }
    for (int i=0; i<otherSparseMatrix.n+1; i++)
        *(colptr+i) = otherSparseMatrix.colptr[i];

    return *this;
}

// unary +
SparseMatrix SparseMatrix::operator+() const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
        new_nzvalues[i] = + this->nzval[i]; 

    // build back cidx
    int* cidx_bis = new int[nnz];
    for (int j = 0; j < this->n; j++) 
    {
        for (int k = colptr[j]; k < colptr[j+1]; k++)
            cidx_bis[k] = j;
    }

    SparseMatrix new_plus_sm(this->nnz, this->rowidx, cidx_bis, new_nzvalues, this->GetSize(1), this->GetSize(2));
    delete[] new_nzvalues; 
    delete[] cidx_bis;
    
    return new_plus_sm;
}

// unary -
SparseMatrix SparseMatrix::operator-() const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
        new_nzvalues[i] = - this->nzval[i];
    
    // build back cidx
    int* cidx_bis = new int[nnz];
    for (int j = 0; j < this->n; j++) 
    {
        for (int k = colptr[j]; k < colptr[j+1]; k++)
            cidx_bis[k] = j;
    }
    SparseMatrix new_minus_sm(this->nnz, this->rowidx, cidx_bis, new_nzvalues, this->GetSize(1), this->GetSize(2));
    delete[] new_nzvalues;
    delete[] cidx_bis;
    
    return new_minus_sm;
}

// scalar multiplication
SparseMatrix SparseMatrix::operator *(double a) const
{
    double *new_nzvalues;
    new_nzvalues = new double[this->nnz];

    for (int i=0; i<this->nnz; i++)
        new_nzvalues[i] = a * this->nzval[i];
    
    // build back cidx
    int* cidx_bis = new int[nnz];
    for (int j = 0; j < this->n; j++) 
    {
        for (int k = colptr[j]; k < colptr[j+1]; k++)
            cidx_bis[k] = j;
    }

    SparseMatrix new_multiply_sm(this->nnz, this->rowidx, cidx_bis, new_nzvalues, this->GetSize(1), this->GetSize(2));
    delete[] new_nzvalues;
    delete[] cidx_bis;
    
    return new_multiply_sm;
}


// product between a vector and a sparse matrix
Vector operator *(const SparseMatrix& m, const Vector& v)
{
    Vector v_product(m.GetSize(1));
    int col = 0;     // keeps track of which column we are in
    int index = 0;   // keeps track inside the given column
    while(index < m.nnz)
    {
        for (int i = 0; index+i < m.colptr[col+1]; i++)
            v_product[m.rowidx[index+i]] += m.nzval[index+i] * v.Read(col);   
        index = m.colptr[col+1];
        col++;
    }
    return v_product;
}