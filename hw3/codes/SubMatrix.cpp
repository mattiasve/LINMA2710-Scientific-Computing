#include "SubMatrix.hpp"
#include <cassert>

double* SubMatrix::ptr(int i, int j)
{
    assert(i >= 1 and i <= numRows);
    assert(j >= 1 and j <= numCols);
    return &(mData[(i-1) * numCols + (j-1)]);
}

SubMatrix::SubMatrix(int numRows, int numCols)
{
    assert(numRows > 0);
    assert(numCols > 0);
    this->numRows = numRows;
    this->numCols = numCols;
    mData = new double[numRows * numCols];
    for (int i = 0; i < numRows * numCols; i ++)
        mData[i] = 0.;
    MPI_Type_vector(numRows, 1, numCols, MPI_DOUBLE, &stride);
    MPI_Type_commit(&stride);
}

SubMatrix::~SubMatrix()
{
    delete[] mData;
    MPI_Type_free(&stride);
}

int SubMatrix::Size(int row) const
{
    assert(row == 1 or row == 2);
    if (row == 1)
        return this->numRows;
    return this->numCols;
}

double SubMatrix::Read(int i, int j) const // 1-based index
{
    assert(i >= 1 and i <= numRows);
    assert(j >= 1 and j <= numCols);
    return mData[(i-1) * numCols + (j-1)];
}

double& SubMatrix::operator()(int i, int j) // 1-based index
{
    assert(i >= 1 and i <= numRows);
    assert(j >= 1 and j <= numCols);
    return mData[(i-1) * numCols + (j-1)];
}

void SubMatrix::SendReceiveRows(int rowsend, int ranksend, int rowrecv, int rankrecv, int tag, MPI_Comm comm){
    MPI_Sendrecv(
        ptr(rowsend, 1), numCols, MPI_DOUBLE, ranksend, tag,  // row send
        ptr(rowrecv, 1), numCols, MPI_DOUBLE, rankrecv, tag,  // row received
        comm, MPI_STATUS_IGNORE
    );
}

void SubMatrix::SendReceiveColumns(int colsend, int ranksend, int colrecv, int rankrecv, int tag, MPI_Comm comm){
    MPI_Sendrecv(
        ptr(1, colsend), 1, stride, ranksend, tag,  // col send 
        ptr(1, colrecv), 1, stride, rankrecv, tag,  // col received
        comm, MPI_STATUS_IGNORE
    );
}
