#include <iostream>
// #include <vector.cpp>

// double inProd = Vector::GetInnerProdOf(a, b);
// Function to allocate memory for a matrix dynamically
double** AllocateMemoryMatrix(int numRows, int numCols)
{
    double** matrix ; // "matrix" is an array of pointers to double, hence "**"
    matrix = new double* [numRows] ;  // allocate memory for an array of pointers (to double) of length numRows
                                      // Actually, the variable matrix contains the address of the pointer matrix[0]
    for (int i=0; i<numRows; i++)
    {
        matrix[i] = new double [numCols] ; // allocate memory for an array of double of length numCols
                                           // Actually, the variable matrix[i] contains the address of matrix[i][0]
    }
    return matrix;
}

// Function to free memory of a matrix
void FreeMatrixMemory(int numRows, double** matrix)
{
    for (int i=0; i<numRows; i++)
    {
        delete [] matrix[i]; // note that numCols is not needed
    }
    delete matrix ;
}

void DisplayMatrix(double** matrix, int numRows, int numCols)
{
    for (int i=0; i<numRows; i++)
    {
        for (int j=0; j<numCols; j++)
        {
            std::cout << matrix[i][j] ;
            if (j==numCols-1)
            std::cout << "\n" ;
        }
    }
    std::cout << "\n" ;   
}

int main()
{
    // Dynamic memory allocation example:
    double** A ;
    A = AllocateMemoryMatrix(5,3); 
    // DisplayMatrix(A, 5, 3) ;
    A[0][1] = 2.0 ;
    A[4][2] = 4.0 ;

    DisplayMatrix(A, 5, 3) ;
    FreeMatrixMemory(5, A) ; // It is crucial to delete dynamically-allocated memory !!
}