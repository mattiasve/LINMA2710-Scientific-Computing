#include "Vector.hpp"
#include "SparseVector.hpp"
#include <cassert>
#include <iostream>

SparseVector::SparseVector() : AbstractVector(0) // empty default constructor
{
    // à priori c'est OK
    nnz = 0; 
    rowidx = NULL; 
    nzval = NULL; 
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
    assert(i < length(*this)); 

    for (int j = 0; j < length(*this); j++)
    {
    if (*(rowidx+j) == i)
        return *(nzval+j);
    }
    return 0.0;
}

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{   
    for(int i=0; i<length(otherVector); i++) // iterate over the size of otherVector
    {
        nzval[i] = otherVector.Read(i) ; // acces value of private member with member function
        // rowidx[i] = otherVector.Read(i) ; // fonctionne pas car Read lis return seulement nzval[i]
    }
    return *this ; 
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    SparseVector SumOfVector; 
    SumOfVector = *this ; 

    return v1; 
}

