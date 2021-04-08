#include <iostream>
#include "myLib.H"

using namespace std;

void myInt::assignA(int a)
{
    a_ = a;
}

void myInt::assignB(int b)
{
    b_ = b;
}

void myInt::sum()
{
    c_ = a_ + b_;
}

void myInt::output()
{
    cout << c_ << endl;
}
