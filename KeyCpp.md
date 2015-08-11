This C++ library provides a MATLAB-like interface to many useful numerical methods. Plotting capability is also provided through a MATLAB-like syntax using Gnuplot.

# Documentation #
The full documentation is listed [here](http://people.tamu.edu/~mons140/). Also provided in the documentation are some basic instructions for compiling and using this library.

# MATLAB/Octave to KeyCpp Syntax Conversion #
Although the goal of this library is to offer a C++ interface similar in syntax to MATLAB/Octave, there are some minor differences. The table below provides equivalent KeyCpp syntax for many common MATLAB/Octave commands. This conversion chart is not complete. If you can't find the conversion you are looking for you can try looking at the [documentation](http://people.tamu.edu/~mons140/) and submitting an issue/feature request if appropriate.

Note: You can omit the `keycpp::` prefix from the following commands by placing `using namespace keycpp;` in the same scope. This shortcut should be used with care as collisions with other libraries are possible.

| **MATLAB/Octave** | **KeyCpp** | Notes |
|:------------------|:-----------|:------|
| `A(1,1)`          | `A(0,0);`  | Indexing starts at 0 in KeyCpp |
| `A(N,N)`          | `A(N-1,N-1);` |       |
| `size(A,1)`       | `keycpp::size(A,1);` |       |
| `size(A,2)`       | `keycpp::size(A,2);` |       |
| `A(:,k)`          | `A.getCol(k-1);` | C++ restricts the use of `:` |
| `A(k,:)`          | `A.getRow(k-1);` |       |
| `A.'`             | `keycpp::transpose(A);` | C++ does not allow overloading `.'` |
| `A'`              | `keycpp::ctranspose(A);` | C++ does not allow overloading `'` |
| `A = zeros(m,n)`  | `keycpp::matrix<double> A = keycpp::zeros(m,n);` | or more simply: `keycpp::matrix<double> A(m,n);` |
| `A = ones(m,n)`   | `keycpp::matrix<double> A = keycpp::ones(m,n);` |       |
| `A.*B`            | `keycpp::times(A,B);` | C++ does not allow overloading `.*` or `./` |
| `A./B`            | `keycpp::rdivide(A,B);` |       |
| `A\b`             | `keycpp::linsolve(A,b);` | `b` is a vector |
| `[V, D] = eig(A,B)` | `std::vector<std::complex<double>> d = keycpp::eig(A,B,&V);` | Non-Hermitian generalized eigenvalue/eigenvector solver uses LAPACK. |
| `x = linspace(0,10,N_x)` | `std::vector<double> x = keycpp::linspace(0.0,10.0,N_x);` |       |
| `x = logspace(1,3,N_x)` | `std::vector<double> x = keycpp::logspace(1.0,3.0,N_x);` | `10 <= x <= 1000` |
| `A = diag([a1, a2, a3])` | `keycpp::matrix<double> A = keycpp::diag({a1, a2, a3});` | `a1`, `a2`, and `a3` are scalar elements of `A` |
| `A = [[a1, a2]; [a3, a4]]` | `keycpp::matrix<double> A = {{a1, a2}, {a3, a4}};` |       |

# Example Code #

## example.cpp ##
```
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
    disp(x);

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
```

**Text Output:**
<br>
1.76667<br>
<br>
-0.333333<br>
<br>
<b>Plot Output</b>
<br>
<br>
<img src='http://people.tamu.edu/~mons140/plot1.png'>
<img src='http://people.tamu.edu/~mons140/plot2.png'>