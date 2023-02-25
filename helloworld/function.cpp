#include <iostream>

// FUNCTIONS
int Add(int a, int b); 
void SetPointedTo2( int * p_i)
{
    *p_i = 2; //stores 2 in memory with adress p_i
}

void SetTo2(int& i)
{
    i = 2;
}

int main()
{
    int z ;
    z = Add(5,3) ;
    std::cout << "The results is " << z << std::endl;

    int c = 1 ;
    int * p_c ;
    p_c = &c ;
    SetPointedTo2(p_c) ;
    std::cout << "c = " << c << "\n" ; // prints : c = 2

    int b = 1;
    SetTo2(b);
    std::cout << "b = " << b << "\n" ;
}

int Add(int a, int b)
{
    int r;
    r = a+b;
    return r;
}