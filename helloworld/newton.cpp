#include <iostream>
#include <cmath>

// implementation of newton iteration
double SolveNewton (double (*pFunc) (double), double (*pFuncPrime) (double), double x)
{
    double step ;
    do 
    {
        step = (*pFunc)(x)/(*pFuncPrime)(x) ; 
        x -= step ;
    } while (fabs(step) > 1e-5) ; 
    return x ; 
}

double f(double x)
{ return x*x - 2 ; }

double fprime(double x)
{return 2*x ; }

// bisection method implementation
double Bissection (double (*pFunc) (double), double a, double b)
{
    double c ;

    while ((b-a)/2. >= 1e-5) // interval is greater than tolerance
    {
        c = (a+b)/2. ; // new mid-point
        if ((*pFunc)(c) == 0.0) // check if mid-point is root
            return c ; 
        
        else if ((*pFunc)(c) * (*pFunc)(a) < 0) // update interval boundaries
            b = c ; 
        else 
            a = c ; 
    }

    return c ; // return final mid-point 
}

int main()
{
    std::cout << "Newton output is : Root = " << SolveNewton(f, fprime, 1) << std::endl; 
    std::cout << "Bissection output is : Root = " << Bissection(f, 0, 2) << std::endl ;
}