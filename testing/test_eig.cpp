#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(eig_test)
{
    double tol = 1e-10;

    keycpp::matrix<double> A = {{1.0, 2.0},
                                {2.0, 1.0}};
                                
    // Expected eigenvalues
    keycpp::vector_k<double> lambda1 = {-1.0, 3.0};
    keycpp::vector_k<double> lambda2 = {3.0, -1.0};
    
    keycpp::matrix<std::complex<double>> vec;
    auto ans = eig(A, &vec);
    
    BOOST_CHECK((abs(keycpp::norm(lambda1-ans)) < tol) || (abs(keycpp::norm(lambda2-ans)) < tol));
    if(abs(keycpp::norm(lambda1-ans)) < tol)
    {
        BOOST_CHECK_CLOSE(real(vec(0,0)/vec(1,0)), -1.0,tol);
        BOOST_CHECK_CLOSE(real(vec(0,1)/vec(1,1)), 1.0,tol);
    }
    else if(abs(keycpp::norm(lambda2-ans)) < tol)
    {
        BOOST_CHECK_CLOSE(real(vec(0,0)/vec(1,0)), 1.0,tol);
        BOOST_CHECK_CLOSE(real(vec(0,1)/vec(1,1)), -1.0,tol);
    }
}

BOOST_AUTO_TEST_SUITE_END()
