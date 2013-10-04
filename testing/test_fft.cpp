#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(fft_test)
{
    // Test fft
    size_t L = 1000;
    double Fs = 1000;
    double T = 1.0/Fs;
    keycpp::vector_k<double> f = Fs*keycpp::linspace(0.0,1.0,L);
    
    keycpp::vector_k<double> t = keycpp::linspace(0.0,(double)(L-1),L)*T;
    keycpp::vector_k<double> y = keycpp::sin(2.0*keycpp::pi*t);
    keycpp::vector_k<std::complex<double>> Y = keycpp::fft(y)*2.0/1000.0;
    
    BOOST_CHECK(abs(abs(Y[1]) - 1.0) < 1e-6);
    BOOST_CHECK(abs(abs(Y[0]) - 0.0) < 1e-6);
    BOOST_CHECK(abs(abs(Y[L-1]) - 1.0) < 1e-6);
    BOOST_CHECK(abs(abs(Y[L-2]) - 0.0) < 1e-6);
    BOOST_CHECK(abs(keycpp::angle(Y[0])) < 1e-6);
    BOOST_CHECK(abs(keycpp::angle(Y[1]) + keycpp::pi/2.0) < 1e-6);
    BOOST_CHECK(abs(keycpp::angle(Y[L-1]) - keycpp::pi/2.0) < 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
