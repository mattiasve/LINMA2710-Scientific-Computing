#include "poisson.hpp"

void jacobi_iteration(const DistributedMatrix &a, const DistributedMatrix &f, DistributedMatrix &b){
    double h;
    int n, i, j;
    n = a.Size(1)-1;
    h = 1/n;
    for (i = a.FirstWriteRow(); i <= a.LastWriteRow(); i++)
    {
        for (j = a.FirstWriteColumn(); j <= a.LastWriteColumn(); j++)
        {
            if (i == 1 or i == n+1 or j == 1 or j == n+1)
                b(i, j) = a.Read(i, j);  // on the boundary
            else
                b(i, j) = (a.Read(i-1, j) + a.Read(i+1, j) + a.Read(i, j-1) + a.Read(i, j+1) - h*h * f.Read(i, j)) / 4; // apply jacobi formula
        }
    }
}

double sum_squares(const DistributedMatrix &a, const DistributedMatrix &b)
{
    int i, j;
    double delta, current_sum, big_sum;
    current_sum = 0.;
    for (i = a.FirstWriteRow(); i <= a.LastWriteRow(); i++)
    {
        for (j = a.FirstWriteColumn(); j <= a.LastWriteColumn(); j++)
        {
            delta = b.Read(i, j) - a.Read(i, j);
            current_sum += delta*delta;
        }
    }
    MPI_Allreduce(&current_sum, &big_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return big_sum;
}