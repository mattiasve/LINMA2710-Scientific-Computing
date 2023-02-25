#include <iostream>

int addition (int a, int b)
{ return (a+b) ; }

int operation (int x, int y, int (*functocall) (int, int))
{
    int g ;
    g = (*functocall)(x, y) ;
    return g ;
}

int main()
{
    int (*p_function) (int a, int b) ; // declare a pointer to a function with two int parameters that returns an int
    p_function = &addition ; // the function pointed to by p_function is addition

    std::cout << (*p_function)(1,2) << std::endl ; 
    std::cout << p_function(1,2) << std::endl ;  // equivalent to previous line 
    std::cout << operation(1, 2, p_function) << std::endl;
    std::cout << operation(1, 2, &addition) << std::endl; // same
    std::cout << operation(1, 2, addition) << std::endl;  // equivalent

}