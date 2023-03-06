#include "SparseVector.hpp"
#include <cassert>
#include <iostream>

SparseVector::SparseVector() : AbstractVector(1) // empty default constructor
{
    nnz = 0; 
    rowidx = NULL; 
    nzval = NULL; 
}

SparseVector::SparseVector(int nnzero, int const *rowindex, double const *nzvalues, int size) : AbstractVector(size)
{
    nnz = nnzero;
    rowidx = new int[nnzero];
    nzval = new double[nnzero]; 
    for(int i=0; i<nnzero; i++)
    {
        *(rowidx+i) = rowindex[i]; 
        *(nzval+i) = nzvalues[i]; 
    }
}

SparseVector::~SparseVector()
{
    delete [] rowidx ;
    delete [] nzval ;
}

double SparseVector::Read(int i) const
{
    assert(i > -1);
    assert(i < length(*this));

    for (int j = 0; j < this->nnz; j++)
    {
    if (*(rowidx+j) == i)
        return *(nzval+j);
    }
    return 0.0;
}

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{   
    this->mSize = otherVector.GetSize(); 
    this->nnz = otherVector.nnz;
    this->nzval = new double[otherVector.nnz]; 
    this->rowidx = new int[otherVector.nnz];
    for (int i=0; i<otherVector.nnz; i++)
    {
        *(nzval+i) = otherVector.nzval[i];
        *(rowidx+i) = otherVector.rowidx[i];
    }
    return *this ; 
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    // create containers that will store values when iterating
    // and will be copied into a correct-size array later
    int *cntn_idx = new int[nnz+v1.nnz]; 
    double *cntn_val = new double [nnz+v1.nnz]; 

    // iterate over non-zero values of both sparse vector
    // k keeps track of the optimal size
    int i=0, j=0, k=0 ; 
    double sum ; 

    while (i < nnz && j<v1.nnz)
    {
        if (*(rowidx+i) == *(v1.rowidx+j))
        // both values are non-zero values
        // compute and store sum
        {
            sum = *(nzval+i) + *(v1.nzval+j);
            cntn_val[k] = sum;
            cntn_idx[k] = *(rowidx+i); 
            k++ ;
            i++ ;
            j++ ; 
        }
        else if (*(rowidx+i) < *(v1.rowidx+j))
        // only one value is a non-zero value
        // only store the relevant value
        {
            cntn_val[k] = *(nzval+i);
            cntn_idx[k] = *(rowidx+i);
            k++ ; 
            i++ ;
        }
        else 
        // likewise, but for the other sparse vetor
        {
            cntn_val[k] = *(v1.nzval+j);
            cntn_idx[k] = *(v1.rowidx+j);
            k++;
            j++;
        }
    }

    while(i < nnz)
    // append the remaining values in the first sparse vector
    {
        cntn_val[k] = *(nzval+i);
        cntn_idx[k] = *(rowidx+i);
        k++ ; 
        i++ ;
    }

    while(j < v1.nnz)
    // append the remaining values in the second vector
    {
        cntn_val[k] = *(v1.nzval+j);
        cntn_idx[k] = *(v1.rowidx+j);
        k++;
        j++;
    }

    // create optimal size arrays and copy data inside it
    int *opti_index = new int[k] ; 
    double *opti_values = new double[k] ;
    std::copy(cntn_idx, cntn_idx+k, opti_index) ; 
    std::copy(cntn_val, cntn_val+k, opti_values) ;

    // create SparseVector containing the sum
    SparseVector sum_sv(k, opti_index, opti_values, length(v1)); 
    return sum_sv; 
}
