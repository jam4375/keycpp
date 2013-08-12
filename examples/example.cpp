#include <iostream>
#include <keycpp/keycpp.h>
using namespace std;
using namespace keycpp;

// Define a class for our ordinary differential equation that we will use later:
class OdeClass
{
    public:
        void operator()(const vector<double> &y,
                        vector<double> &dy,
                        const double t)
        {
            dy[0] = y[1]*y[2];
            dy[1] = -y[0]*y[2];
            dy[2] = -0.51*y[0]*y[1];
        }
};

int main(int argc, char** argv)
{
    // First, lets create some data: y1 = t^2 and y2 = t^3
    vector<double> t = linspace(-2.0,2.0,100);
    vector<double> y1 = times(t,t);
    vector<double> y2 = times(t,times(t,t));

    // Now, lets plot the data we just created:
    Figure h1;
    h1.plot(t,y1,"b-","linewidth",2);
    h1.hold_on();
    h1.plot(t,y2,"r--","linewidth",2);
    h1.grid_on();
    h1.xlabel("t");
    h1.ylabel("y");
    h1.legend({"y1 = t^2","y2 = t^3"});
    set(h1,"fontsize",14);

    // This is how to solve linear equations of the form Ax = b:
    matrix<double> A = {{1.0, 2.0},
                        {1.0,-1.0}};
    vector<double> b = {1.1,
                        2.1};
    vector<double> x = linsolve(A,b);
    // Print the result to the screen:
    cout << x;

    // Now lets do something a little more complicated, solve an
    // ordinary differential equation (ODE):
    // y(1)' = y(2)*y(3);
    // y(2)' = -y(1)*y(3);
    // y(3)' = 0.51*y(1)*y(2);
    // With initial conditions at t = 0: y(1) = 0; y(2) = 1; y(3) = 1;
    OdeClass myOde;
    vector<double> t2 = linspace(0.0,12.0,100);
    vector<double> ICs = {0.0, 1.0, 1.0};
    matrix<double> y = ode45(myOde, t2, ICs);
    
    // Now that we have solved the ODE, lets plot the results:
    Figure h2;
    h2.plot(t2,y.getCol(0),"-");
    h2.hold_on();
    h2.plot(t2,y.getCol(1),"-.");
    h2.plot(t2,y.getCol(2),"x");
    h2.xlabel("t");
    h2.ylabel("y");
    set(h2,"fontsize",14);
    h2.title("ODE Solution");

    return 0;
}
