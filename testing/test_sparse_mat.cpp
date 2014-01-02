#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(sparse_matrix_test)
{
    double tol = 1e-10;
    
    // Test indexing
    {
        keycpp::matrix<double,2,1> A1 = {{1.0, 2.0},
                                         {3.0, 1.0}};
        BOOST_CHECK(abs(A1(0,0) - 1.0) < tol);
        BOOST_CHECK(abs(A1(1,0) - 3.0) < tol);
        BOOST_CHECK(abs(A1(0,1) - 2.0) < tol);
        BOOST_CHECK(abs(A1(1,1) - 1.0) < tol);
        
        BOOST_CHECK(abs(A1(0) - 1.0) < tol);
        BOOST_CHECK(abs(A1(1) - 3.0) < tol);
        BOOST_CHECK(abs(A1(2) - 2.0) < tol);
        BOOST_CHECK(abs(A1(3) - 1.0) < tol);
        
        // Now modify some things:
        A1(0,0) = 5.5;
        A1(1,0) = 4.4;
        BOOST_CHECK(abs(A1(0,0) - 5.5) < tol);
        BOOST_CHECK(abs(A1(1,0) - 4.5) < tol);
        BOOST_CHECK(abs(A1(0,1) - 2.0) < tol);
        BOOST_CHECK(abs(A1(1,1) - 1.0) < tol);
        
        A1(3) = 5.5;
        A1(2) = 4.4;
        
        BOOST_CHECK(abs(A1(0) - 5.5) < tol);
        BOOST_CHECK(abs(A1(1) - 4.5) < tol);
        BOOST_CHECK(abs(A1(2) - 4.5) < tol);
        BOOST_CHECK(abs(A1(3) - 5.5) < tol);
        
        keycpp::matrix<double,2,1> A2 = {{1.0, 2.0},
                                         {3.0, 1.0}};
        keycpp::matrix<double,2,1> B2 = {3.0, 1.0};
        keycpp::matrix<double,2,1> C2 = {2.0, 1.0};
        
        BOOST_CHECK(abs(keycpp::norm(B2 - A2.row(1))) < tol);
        BOOST_CHECK(abs(keycpp::norm(C2 - keycpp::transpose(A2.col(1)))) < tol);
        
        keycpp::matrix<double,2,1> A3 = {{1.0, 2.0},
                                         {3.0, 1.0}};
        keycpp::matrix<double,2,1> B3 = {5.0, 2.0};
        keycpp::matrix<double,2,1> C3 = {4.0, 3.0};
        
        A3.row(1) = B3;
        BOOST_CHECK(abs(keycpp::norm(B3 - A3.row(1))) < tol);
        A3.col(1) = C3;
        BOOST_CHECK(abs(keycpp::norm(C3 - keycpp::transpose(A3.col(1)))) < tol);
    }
        
    // Test span
    {
        keycpp::matrix<double,2,1> A1 = {{1.0, 2.0},
                                         {3.0, 1.0}};
        keycpp::matrix<double,2,1> B1 = {1.0, 2.0};
        keycpp::matrix<double,2,1> C1 = {1.0, 3.0};
        BOOST_CHECK(abs(keycpp::norm(A1 - A1(keycpp::span(),keycpp::span()))) < tol);
        BOOST_CHECK(abs(keycpp::norm(A1 - A1(keycpp::span(0,1),keycpp::span(0,1)))) < tol);
        BOOST_CHECK(abs(keycpp::norm(A1 - A1(keycpp::span(0,1,1),keycpp::span(0,1,1)))) < tol);
        
        BOOST_CHECK(abs(keycpp::norm(B1 - A1(0,keycpp::span(0,1,1)))) < tol);
        BOOST_CHECK(abs(keycpp::norm(B1 - A1(keycpp::span(0,1,0),keycpp::span(0,1,1)))) < tol);
        
        BOOST_CHECK(abs(keycpp::norm(keycpp::transpose(C1) - A1(keycpp::span(0,1,1),0))) < tol);
        BOOST_CHECK(abs(keycpp::norm(keycpp::transpose(C1) - A1(keycpp::span(0,1,1),keycpp::span(0,1,0)))) < tol);
    }

    // Test matrix double-double multiplication
    {
        keycpp::matrix<double,2,1> A1 = {{1.0, 2.0},
                                         {3.0, 1.0}};
        keycpp::matrix<double,2,1> B1 = {{1.0, -2.0},
                                         {-2.0, 4.0}};
        keycpp::matrix<double,2,1> C1 = {{-3.0, 6.0},
                                         {1.0, -2.0}};
        
        BOOST_CHECK(abs(keycpp::norm(C1 - A1*B1)) < tol);
        BOOST_CHECK(abs(keycpp::norm(C1.row(0) - A1.row(0)*B1)) < tol);
        BOOST_CHECK(abs(keycpp::norm(C1.row(0).col(0) - A1.row(0)*B1.col(0))) < tol);
        
        keycpp::matrix<double,2,1> A2 = {1.0, 2.0};
        keycpp::matrix<double,2,1> B2 = {{1.0, -2.0},
                                         {-2.0, 4.0}};
        keycpp::matrix<double,2,1> C2 = {-3.0, 6.0};
        
        BOOST_CHECK(abs(keycpp::norm(C2 - A2*B2)) < tol);
    }
    
    // Test matrix complex<double>-double multiplication
    {
        std::complex<double> i = std::complex<double>(0.0,1.0);
        keycpp::matrix<std::complex<double>,2,1> A1 = {{i*1.0, i*2.0},
                                                       {i*3.0, i*1.0}};
        keycpp::matrix<double,2,1> B1 = {{1.0, -2.0},
                                         {-2.0, 4.0}};
        keycpp::matrix<std::complex<double>,2,1> C1 = {{-i*3.0, i*6.0},
                                                       {i*1.0, -i*2.0}};
        
        BOOST_CHECK(abs(keycpp::norm(C1 - A1*B1)) < tol);
        BOOST_CHECK(abs(keycpp::norm(C1.row(0) - A1.row(0)*B1)) < tol);
        
        keycpp::matrix<std::complex<double>,2,1> A2 = {i*1.0, i*2.0};
        keycpp::matrix<double,2,1> B2 = {{1.0, -2.0},
                                         {-2.0, 4.0}};
        keycpp::matrix<std::complex<double>,2,1> C2 = {-i*3.0, i*6.0};
        
        BOOST_CHECK(abs(keycpp::norm(C2 - A2*B2)) < tol);
    }
    
    // Test matrix complex<double>-complex<double> multiplication
    {
        std::complex<double> i = std::complex<double>(0.0,1.0);
        keycpp::matrix<std::complex<double>,2,1> A1 = {{i*1.0, i*2.0},
                                                       {i*3.0, i*1.0}};
        keycpp::matrix<std::complex<double>,2,1> B1 = {{i*1.0, -i*2.0},
                                                       {-i*2.0, i*4.0}};
        keycpp::matrix<std::complex<double>,2,1> C1 = {{-i*i*3.0, i*i*6.0},
                                                       {i*i*1.0, -i*i*2.0}};
        
        BOOST_CHECK(abs(keycpp::norm(C1 - A1*B1)) < tol);
        BOOST_CHECK(abs(keycpp::norm(C1.row(0) - A1.row(0)*B1)) < tol);
        BOOST_CHECK(abs(keycpp::norm(C1.row(0).col(0) - A1.row(0)*B1.col(0))) < tol);
        
        keycpp::matrix<std::complex<double>,2,1> A2 = {i*1.0, i*2.0};
        keycpp::matrix<std::complex<double>,2,1> B2 = {{i*1.0, -i*2.0},
                                                       {-i*2.0, i*4.0}};
        keycpp::matrix<std::complex<double>,2,1> C2 = {-i*i*3.0, i*i*6.0};
        
        BOOST_CHECK(abs(keycpp::norm(C2 - A2*B2)) < tol);
    }
}

BOOST_AUTO_TEST_SUITE_END()
