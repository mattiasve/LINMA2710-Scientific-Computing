#include "Vector.hpp"
#include "SparseVector.hpp"
#include <cassert>

SparseVector::SparseVector(){}  // empty default constructor

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

// SparseVector& operator=(const SparseVector& otherVector)
// {
    
// }

int main()
{
    return 0;
}