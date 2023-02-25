#include <iostream>
#include <cassert> // for the assert function

// overloading operator + and operator =
class CVector //define class vector
{
public:
    int x, y ; // declare two member variables
    CVector() {}; // default constructor
    CVector(int a, int b) : x(a), y(b) {} ; // declare and define specialized constructor using member initialization
    CVector operator+ (const CVector&) const ;
    CVector& operator= (const CVector& param) ;
}; // end of class definition

CVector CVector::operator+ (const CVector& param) const
/* define member function "operator+", to which an object (locally labelled "param") of
class CVector is passed by reference, and which returns an object of class CVector */
{
CVector temp ; // create object "temp" of class CVector using default constructor
temp.x = x + param.x ; 
temp.y = y + param.y ;
return temp ; 
}

CVector& CVector::operator=(const CVector& param)
{
    x = param.x ;
    y = param.y ;
    return *this; // * is the dereference operator !!
}

/* -------------------------------------------------------------------------------------- */

// Overlaoding () operator : 1-based indexing
class Vector
{
private: 
    double* mData ; 
    int mSize ;

public:
    Vector(int size) ; // specialized constructor
    ~Vector(); // destructor: needed for the delete that matches the new in the def of constructor
    int GetSize() const ;
    double Read(int i) const ; 
    double& operator() (int i); // overloading operator()
    friend double operator*(Vector& a, Vector& b) ; 
    friend int length(const Vector& v); // friend function "length"
    /* In C++, a friend function is a function that is not a member of a class, but is allowed to 
    access its private and protected members. In this case, length is defined outside of the Vector 
    class, but has access to its private data members because it is declared as a friend. */
    double GetInnerProdOf(const Vector& v, const Vector& w);
}; // end of class definition

Vector::Vector(int size)
{
    mData = new double[size] ; // memory allocation od an array of double of size "size"
    mSize = size ; 
}

Vector::~Vector() 
{
    delete [] mData ; 
}

double Vector::Read(int i) const {return mData[i] ;}

int Vector::GetSize() const {return mSize;}

double& Vector::operator()(int i)
{
    assert(i >= 1); // assert will terminate the program if its argument is not satisfied
    assert(i <= mSize); 
    return mData[i-1]; 
}

int length(const Vector& v) {return v.mSize;} // def of friend function length

double Vector::GetInnerProdOf(const Vector& v, const Vector& w)
{
    double sum ; 
    sum = 0.0 ;

    for (int i= 0; i<v.GetSize(); i++)
    {
        sum = sum + v.Read(i) * w.Read(i);
    }
    return sum ; 
}

double operator* (Vector& a, Vector& b) 
{
    return a.GetInnerProdOf(a, b) ;
}

int main()
{
    CVector u (3,1) ; 
    CVector v (1,2) ; 
    CVector result1, result2, r ; 
    result1 = u.operator+(v) ;// oprerator+ is just like a normal member function of CVector
    result2 = u+v ;// but is can be called with a simplier syntax
    r.operator=(u) ;
    r = u ; // equivalent to previous line
    std::cout << result1.x << "," << result1.y << "\n";
    std::cout << result2.x << "," << result2.y << "\n";
    std::cout << u.x << "," << u.y << "\n";
    std::cout << "u.x, u.y = " << u.x << "," << u.y << "\n" ;
    std::cout << "r.x, r.y = " << r.x << "," << r.y << "\n\n" ;

    Vector w(2), a(2) ; 
    // Vector v{2} : C++11, equivalent as above
    std::cout << "w(1), w(2) = " << w(1) << "," << w(2) << "\n" ; 
    w(1) = 3 ; 
    w(2) = 4 ; 
    a(1) = 2 ;
    a(2) = 3 ;
    std::cout << "w(1), w(2) = " << w(1) << "," << w(2) << "\n" ; 
    std::cout << "a(1), w(2) = " << a(1) << "," << a(2) << "\n" ; 
    std::cout << "w.GetSize() = " << w.GetSize() << "\n" ;
    std::cout << "length(w) = " << length(w) << "\n" ; 
    std::cout << "inner product of w.a = " << w.GetInnerProdOf(w, a) << "\n" ; 

    return 0 ;
}