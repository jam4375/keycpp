#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(sprintf_test)
{
    // Test sprintf
    double a = 2.223;
    BOOST_CHECK(keycpp::sprintf("a = %0.2f",a).compare("a = 2.22") == 0);
}

BOOST_AUTO_TEST_SUITE_END()
