#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(fft_test)
{
    // Test fft
    {
        size_t L = 1000;
        double Fs = 1000;
        double T = 1.0/Fs;
        keycpp::matrix<double,2> f = Fs*keycpp::linspace(0.0,1.0,L);
        
        keycpp::matrix<double,2> t = keycpp::linspace(0.0,(double)(L-1),L)*T;
        keycpp::matrix<double,2> y = keycpp::sin(2.0*keycpp::pi*t);
        keycpp::matrix<std::complex<double>,2> Y = keycpp::fft(y)*2.0/1000.0;
        
        BOOST_CHECK(abs(abs(Y(1)) - 1.0) < 1e-6);
        BOOST_CHECK(abs(abs(Y(0)) - 0.0) < 1e-6);
        BOOST_CHECK(abs(abs(Y(L-1)) - 1.0) < 1e-6);
        BOOST_CHECK(abs(abs(Y(L-2)) - 0.0) < 1e-6);
        BOOST_CHECK(abs(keycpp::angle(Y(0))) < 1e-6);
        BOOST_CHECK(abs(keycpp::angle(Y(1)) + keycpp::pi/2.0) < 1e-6);
        BOOST_CHECK(abs(keycpp::angle(Y(L-1)) - keycpp::pi/2.0) < 1e-6);
    }
    
    // Test ifft
    {
        size_t L = 1000;
        double Fs = 1000;
        double T = 1.0/Fs;
        keycpp::matrix<double,2> f = Fs*keycpp::linspace(0.0,1.0,L);
        
        keycpp::matrix<double,2> t = keycpp::linspace(0.0,(double)(L-1),L)*T;
        keycpp::matrix<double,2> y = keycpp::sin(2.0*keycpp::pi*t);
        keycpp::matrix<std::complex<double>,2> Y = keycpp::ifft(keycpp::fft(y));
        
        BOOST_CHECK(keycpp::norm(keycpp::abs(Y-y)) < 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()
