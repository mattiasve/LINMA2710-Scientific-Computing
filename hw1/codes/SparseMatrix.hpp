#ifndef SparseMatrixHEADERDEF
#define SparseMatrixHEADERDEF
#include "SparseVector.hpp"
#include "Vector.hpp"
#include "SparseVector.hpp"
class SparseMatrix   // ELL format Sparse Matrix
{
    private:
        int m;          // Number of rows
        int n;          // Number of columns
        int nnz;        // number of non-zero values
        int *colptr;    // points to the starting indices of the columns in rowidx and nzval
        int *rowidx;    // row indices
        double * nzval; // Stored values
    
    public:
        SparseMatrix(int numRows, int numCols);
        SparseMatrix(int m, int n, int nnz, int* colptr, int* rowidx, double* nzval);
        SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int m, int n);
        ~SparseMatrix();
        int GetSize(int) const;
        //overloaded assignment operator
        SparseMatrix& operator=(const SparseMatrix& otherSparseMatrix);
        SparseMatrix operator+() const; // unary +
        SparseMatrix operator-() const; // unary -
        
        // scalar multiplication
        SparseMatrix operator *(double a) const;
        // declare vector multiplication friendship
        friend Vector operator *(const SparseMatrix& m,
                const Vector& v);
};
// prototype signatures for friend operators
Vector operator *(const SparseMatrix& m, const Vector& v);
#endif
