#include "Vector.hpp"
#include "SparseVector.hpp"
#include <cassert>

//SparseVector::SparseVector(){}  // empty default constructor

SparseVector::SparseVector(int nnz, int const *rowidx, double const *nzval, int size)
{
    nzVector = new double[nnz] ;

    for (int i=0; i<nnz; i++)
    {
        nzVector[rowidx[i]] = nzval[i] ; 
    }
}

SparseVector::~SparseVector()
{
    delete [] nzVector ;
}

double SparseVector::Read(int i) const
{
    assert(i > -1);
    assert(i < GetSize());
    return nzVector[i]; 
}

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{
    SparseVector copyObj; 
    copyObj = otherVector; 
    return copyObj; 
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    assert(GetSize() == v1.GetSize()); 
    SparseVector v0(); // comment construire le vecteur ?? 
    for (int i=0; i<GetSize(); i++)
    {
        v0[i] = nzVector[i] + v1.nzVector[i]; 
    }
    return v0; 
}

int main()
{
    return 0;
}