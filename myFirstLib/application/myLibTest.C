#include <istream>
#include "myLib.H"

using namespace std;

int main()
{
    myInt myClass;

    myClass.assignA(4);

    myClass.assignB(10);

    myClass.output();

    myClass.sum();

    myClass.output();

    return 0;

}
