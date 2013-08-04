#include <iostream>
#include "keycpp.h"
using namespace std;
using namespace keycpp;


int main(int argc, char** argv)
{
	vector<double> x = linspace(-2.0,2.0,100);
	vector<double> y1 = etimes(x,x);
	vector<double> y2 = etimes(x,etimes(x,x));
	Figure h;
	h.plot(x,y1);
	h.hold_on();
	h.plot(x,y2);
	h.grid_on();
	h.xlabel("x");
	h.ylabel("y");
	h.legend({"y1 = x^2","y2 = x^3"});
	
	return 0;
}
