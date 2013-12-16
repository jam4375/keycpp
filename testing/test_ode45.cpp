#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

// Define a class for our ordinary differential equation that we will use later:
class OdeClass
{
    public:
        void operator()(const keycpp::matrix<double,2> &y,
                        keycpp::matrix<double,2> &dy,
                        const double t)
        {
            dy(0) = t;
            dy(0) = y(1);
            dy(1) = -y(0);
        }
};

BOOST_AUTO_TEST_CASE(ode45_test)
{
    OdeClass myOde;
    keycpp::matrix<double,2> t = keycpp::linspace(0.0,12.0,100);
    keycpp::matrix<double,2> ICs = {1.0, 0.0};
    keycpp::matrix<double> y = keycpp::ode45(myOde, t, ICs);
    
    for(size_t ii = 0; ii < t.size(1); ii++)
    {
        BOOST_CHECK_CLOSE(y(ii,0),cos(t(ii)),0.01);
    }
}

BOOST_AUTO_TEST_SUITE_END()
