#include <iostream> // for cout
#include <math.h> // for pow
using namespace std;

double f(double x) { // the function that we want to integrate
    return pow(x,2); // returns x^2
}

int main() {
    double a = 0; // left endpoint of the integration interval
    double b = 1; // right endpoint
    int n = 1e3; // number of subintervals
    double s = 0;
    double h = (b-a)/n; // size of the subintervals
    double x_minus, x_plus, x;
    for (int i=1; i<=n; i++)
    {
        x_minus = a + (i-1)*h;
        x_plus = a + i*h;
        x = (x_minus + x_plus)/2;
        s = s + 1./6.*f(x_minus) + 4./6.*f(x) + 1./6.*f(x_plus);
    }
    s = h*s;
    cout << s << endl; // display result
}
