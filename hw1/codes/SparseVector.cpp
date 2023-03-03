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
    this->nzval = otherVector.nzval; 
    this->rowidx = otherVector.rowidx;

    return *this ; 
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    for (int i=0; i< std::max(this->nnz, v1.nnz); i++) // iterate over maximum value between this.nnz and v1.nnz
    {
        *this->nzval = this->Read(i) + v1.Read(i); 
    }

    return *this; 
}

// int main()
// {
//     return 0; 
// }