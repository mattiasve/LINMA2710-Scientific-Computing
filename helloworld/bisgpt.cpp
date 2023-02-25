#include <iostream>
#include <cmath>

using namespace std;

// Function to find the root using the bisection method
double bisection(double a, double b, double eps, double (*f)(double))
{
    double c;
    while ((b-a) >= eps)
    {
        // Find midpoint
        c = (a+b)/2;

        // Check if midpoint is root
        if (f(c) == 0.0)
            return c;

        // Decide which half to repeat
        else if (f(c)*f(a) < 0)
            b = c;
        else
            a = c;
    }

    // Return midpoint of the final interval
    return (a+b)/2;
}

// Example function to find the root of
double example_function(double x)
{
    return x*x - 2;
}

int main()
{
    double a = 0, b = 2, eps = 0.0001;
    double root = bisection(a, b, eps, &example_function);
    cout << "Root: " << root << endl;
    return 0;
}
