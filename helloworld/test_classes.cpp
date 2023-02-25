#include <iostream>
using namespace std;
// using namespace std;
class Rectangle // Rectangle is the identifier of the class 
{ // begin body of class declaration
    public: // begin public section
        Rectangle(); // declare default constructor 
        Rectangle(int,int); // declare specialized constructor
        void SetWidth(int);
        int GetArea(void) {return (mWidth*mHeight);}
        double id;
        double getId(void);
    private: // Priveta section
        int mWidth, mHeight;
};
Rectangle::Rectangle() //define default constructor
{
    mWidth = 5;
    mHeight = 5;
    id = 2.0;
}
Rectangle::Rectangle(int a, int b) // define specialized constructor 
{
    mWidth = a;
    mHeight = b; 
}

void Rectangle::SetWidth(int a) // define member function SetWidth 
{
    mWidth = a; 
}
double Rectangle::getId(void) {
   return id*2.0;
}
int main()
{
    Rectangle rect (3,4);  //Using specialized constructor
    Rectangle rectb;       //Using default constructor
    Rectangle *p_rect ;     // declare pointer p_rect to an object of class Rectangle
    p_rect = new Rectangle(5,6) ; // memory allocation and initialization of object pointed to by p_rect
    rectb.SetWidth(4);    
    std::cout<< "ID_=" << rectb.getId() << std::endl;
    std::cout << "rect_area:_" << rect.GetArea() << std::endl;
    std::cout << "rectb␣area:␣" << rectb.GetArea() << std::endl;
    delete p_rect ; 
    return 0;

}