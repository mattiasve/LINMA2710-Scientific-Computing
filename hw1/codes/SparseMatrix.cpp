#include "SparseMatrix.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>
#include <vector>

// zero sparse matrix contructor
SparseMatrix::SparseMatrix(int numRows, int numCols)
{
    m = numRows; 
    n = numCols; 
    nnz = 0; 
    
    colptr = new int[numRows+1];
    for (int i=0; i<numRows+1; i++)
        colptr[i] = 0;

    rowidx = new int[nnz];
    nzval = new double[nnz];

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

    // copy data inside arrays
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
    this->m = size1;
    this->n = size2;
    this->nnz = nnz;
    this->nzval = new double[nnz];
    this->rowidx = new int[nnz];
    this->colptr = new int[size2+1];

    /*
    In this following part I tried to implement a method to sort cidx in a way that 
    still matches with ridx and nzval. The methods works but does not yield the right sparse matrix.
    Values in the same column are switched 
    // ===== SORTING =====
    struct Triple 
    {
      int ri;
      int ci;
      double nzv;
    };

    // Create an array of struct that combines ri, ci, and nzv values
    Triple triples[nnz];
    for (int i = 0; i < nnz; i++) 
    {
        triples[i].ri  = *(ridx+i);
        triples[i].ci  = *(cidx+i);
        triples[i].nzv = *(nzval+i);
    }

    // Sort the array of struct based on the ci value
    std::sort(triples, triples + 7,
                [](const Triple& lhs, const Triple& rhs) 
                { return lhs.ci < rhs.ci; });

    // Extract the sorted ri, ci, and nzv arrays from the sorted array of struct
    int sorted_ridx[nnz];
    int sorted_cidx[nnz];
    double sorted_nzval[nnz];
    for (int i = 0; i < nnz; i++) 
    {
        sorted_ridx[i]  = triples[i].ri  ;
        sorted_cidx[i]  = triples[i].ci  ;
        sorted_nzval[i] = triples[i].nzv ;
    }

    // Print the sorted arrays for verification
    for (int i = 0; i < 7; ++i) {
    std::cout << "ridx[" << i << "] = " << sorted_ridx[i] << ", cidx[" << i << "] = " << sorted_cidx[i] << ", nzval[" << i << "] = " << sorted_nzval[i] << std::endl;
    }
    // ===== END SORTING =====
    */

    // fill colptr with zeros
    for (int i=0; i<size2+1; i++)
        this->colptr[i] = 0; 

    // increment colptr by 1 if it corresponds to an column index of an non-zero value
    for (int i=0; i<nnz; i++)
        this->colptr[cidx[i]+1]++;
    
    // compute cumulative sum of non-zero value
    // allows us to dertermine the position in the sparse matrix
    for (int i=1; i<size2+1; i++)
        this->colptr[i] += this->colptr[i-1];

    // copy relative elements in the right position
    for (int i=0; i<nnz; i++)
    {
        int index = this->colptr[cidx[i]]; 
        this->nzval[index]  = nzval[i];
        this->rowidx[index] = ridx[i];
        this->colptr[cidx[i]]++; // increment column pointer for the next non-zero value in the same column
    }

    // reset column pointer to the start of each column
    for (int i=size2; i>0; i--)
        this->colptr[i] = this->colptr[i-1];
    this->colptr[0] = 0;
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
    this->colptr = new int[otherSparseMatrix.n+1];
    this->rowidx = new int[otherSparseMatrix.n];
    this->nzval = new double[otherSparseMatrix.n];

    for (int i=0; i<otherSparseMatrix.n; i++)
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
    return new_multiply_sm;
}

// product between sparse matrix and vector
Vector operator *(const SparseMatrix& m, const Vector& v)
{
    assert(v.GetSize() == m.GetSize(2));

    // create vector to store solution
    Vector sol_v(m.GetSize(1));

    for (int j=0; j<m.GetSize(2); j++)
    {
        // get limits of column j
        int start = m.colptr[j];
        int end = m.colptr[j+1];

        // iterate over the non-zero values in column j
        for (int k=start; k<end; k++)
        {
            //std::cout << " yooo " << std::endl;
            int i = m.rowidx[k];
            double value = m.nzval[k];

            // std::cout << "nzval[j] = " <<  m.nzval[k] << std::endl;
            // std::cout << "v.Read(j) = " << v.Read(j) << std::endl;

            sol_v[i] += value * v.Read(j);
        }
    }
    return sol_v;
}