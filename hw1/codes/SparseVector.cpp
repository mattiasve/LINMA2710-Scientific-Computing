#include "Vector.hpp"
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
    nnz = otherVector.nnz;
    nzval = otherVector.nzval; 
    rowidx = otherVector.rowidx;

    // std::cout << "this->nnz = " << this->nnz << std::endl;
    // for (int i=0; i<nnz; i++)
    //     std::cout << nzval[i] << std::endl;
    // for (int i=0; i<nnz; i++)
    //     std::cout << rowidx[i] << std::endl;

    return *this ; 
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    // SparseVector SumOfVector; 
    // SumOfVector = *this ; 

    return v1; 
}

