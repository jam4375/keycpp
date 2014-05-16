#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(eigs_test)
{
    double tol = 1e-13;
    
    {
        keycpp::matrix<std::complex<double>> Evecs, Evals;
        std::complex<double> i = std::complex<double>(0.0,1.0);
        
        int n = 100;
        keycpp::matrix<std::complex<double>> A = keycpp::rand(n,n) + i*keycpp::rand(n,n);

        int nev = 3; // The number of values to calculate

        Evals = keycpp::eigs(A, nev, "LM", &Evecs);
        BOOST_CHECK((abs(keycpp::norm(A*Evecs - Evecs*Evals)) < tol));
    }
    
    {
        keycpp::matrix<std::complex<double>> Evecs, Evals;
        std::complex<double> i = std::complex<double>(0.0,1.0);
        
        keycpp::matrix<std::complex<double>> A = {{1.0*i, 2.0,1.1},
                                     {0.5*i, 0.5, -0.5},
                                     {-0.3, 4.0,0.2}};
        keycpp::matrix<std::complex<double>> B = {{3.0, -1.0,10.0},
                                     {-1.0, 2.0, -1.0},
                                     {1.0, -1.0,2.0}};

        int nev = 1; // The number of values to calculate

        Evals = keycpp::eigs(A, B, nev, "LM", &Evecs);
        
        BOOST_CHECK((abs(keycpp::norm(A*Evecs - B*Evecs*Evals)) < tol));
    }
    
    /*// Test for sparse matrix
    {
        keycpp::matrix<std::complex<double>> Evecs, Evals;
        std::complex<double> i = std::complex<double>(0.0,1.0);
        
        keycpp::matrix<std::complex<double>> A = {{1.0*i, 2.0,1.1},
                                     {0.5*i, 0.5, -0.5},
                                     {-0.3, 4.0,0.2}};
                                     
        keycpp::matrix<std::complex<double>,2,SPARSE_MATRIX> B = sparse(A);

        int nev = 1; // The number of values to calculate

        Evals = keycpp::eigs(B, nev, "LM", &Evecs);
        
        BOOST_CHECK((abs(keycpp::norm(B*Evecs - Evecs*Evals)) < tol));
    }*/
}

BOOST_AUTO_TEST_SUITE_END()
