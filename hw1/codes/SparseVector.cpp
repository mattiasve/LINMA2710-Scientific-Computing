#include "Vector.hpp"
#include "SparseVector.hpp"
#include <cassert>
#include <iostream>

SparseVector::SparseVector() : AbstractVector(0) // empty default constructor
{
    nnz = 0; 
}

SparseVector::SparseVector(int nnzero, int const *rowindex, double const *nzvalues, int size) : AbstractVector(size)
{
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
    assert(i < GetSize());
    return *(nzval+i);  
}


SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{
    rowidx = otherVector.rowidx ; 
    nzval = otherVector.nzval ;
    return *this ; 
}

int main()
{
    return 0;
}