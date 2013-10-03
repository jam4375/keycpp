#define BOOST_TEST_MODULE "C++ Unit Tests for KeyCpp"
#include <boost/test/included/unit_test.hpp>

/*#define BOOST_TEST_MODULE vector_k_test
#define BOOST_TEST_MODULE keycpp_test
#include <boost/test/included/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_CASE(vector_k_test)
{
    {
        keycpp::vector_k<int> v1;
        BOOST_CHECK_EQUAL(v1.capacity(),0);
        BOOST_CHECK_EQUAL(v1.size(),0);
        BOOST_CHECK(v1.empty());
        BOOST_CHECK(v1.begin() == nullptr);
        BOOST_CHECK(v1.end() == nullptr);
        
        keycpp::vector_k<int> v2(10);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        int* it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 5678;
        v2[9] = 8765;
        BOOST_CHECK_EQUAL(*v2.begin(),5678);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),8765);
        
        v2.resize(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),20);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        const int* it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        const int* it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
        
        v2.push_back(3333);
        BOOST_CHECK(v2.capacity() >= 21);
        BOOST_CHECK_EQUAL(v2.size(),21);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),3333);
        it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),3333);
        
        v2.pop_back();
        BOOST_CHECK(v2.capacity() >= 21);
        BOOST_CHECK_EQUAL(v2.size(),20);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        
        const int aa = v2[0];
        const int bb = v2[19];
        BOOST_CHECK_EQUAL(aa,1111);
        BOOST_CHECK_EQUAL(bb,2222);
        
        keycpp::vector_k<int> v3(5);
        keycpp::vector_k<int> v4;
        int ii = 0;
        for(keycpp::vector_k<int>::iterator it = v3.begin(); it != v3.end(); ++it)
        {
            *it = ii;
            ii++;
        }
        for(keycpp::vector_k<int>::iterator it = v3.begin(); it != v3.end(); ++it)
        {
            v4.push_back(*it);
        }
        BOOST_CHECK_EQUAL(v4.size(),v3.size());
        for(keycpp::vector_k<int>::size_type jj = 0; jj < v4.size(); jj++)
        {
            BOOST_CHECK_EQUAL(v4[jj],v3[jj]);
        }
        
        std::vector<int> v5(50);
        keycpp::vector_k<int> v6;
        ii = 0;
        for(std::vector<int>::iterator it = v5.begin(); it != v5.end(); ++it)
        {
            *it = ii;
            ii++;
        }
        for(std::vector<int>::iterator it = v5.begin(); it != v5.end(); ++it)
        {
            v6.push_back(*it);
        }
        BOOST_CHECK_EQUAL(v5.size(),v6.size());
        for(keycpp::vector_k<int>::size_type jj = 0; jj < v6.size(); jj++)
        {
            BOOST_CHECK_EQUAL(v6[jj],v5[jj]);
        }
        
        keycpp::vector_k<int> v7(50,27);
        keycpp::vector_k<int> v8;
        for(keycpp::vector_k<int>::iterator it = v7.begin(); it != v7.end(); ++it)
        {
            BOOST_CHECK_EQUAL(*it,27);
            v8.push_back(*it);
        }
        BOOST_CHECK_EQUAL(v8.size(),v7.size());
        for(keycpp::vector_k<int>::size_type jj = 0; jj < v8.size(); jj++)
        {
            BOOST_CHECK_EQUAL(v8[jj],v7[jj]);
        }
        
        keycpp::vector_k<int> v9 = {27, 27, 27, 27};
        BOOST_CHECK_EQUAL(v9.size(),4);
        keycpp::vector_k<int> v10(v9);
        for(keycpp::vector_k<int>::iterator it = v9.begin(); it != v9.end(); ++it)
        {
            BOOST_CHECK_EQUAL(*it,27);
        }
        BOOST_CHECK_EQUAL(v10.size(),v9.size());
        for(keycpp::vector_k<int>::size_type jj = 0; jj < v10.size(); jj++)
        {
            BOOST_CHECK_EQUAL(v10[jj],v9[jj]);
        }
        
        std::vector<int> v11(50,27);
        keycpp::vector_k<int> v12(v11);
        BOOST_CHECK_EQUAL(v11.size(),v12.size());
        for(keycpp::vector_k<int>::size_type jj = 0; jj < v12.size(); jj++)
        {
            BOOST_CHECK_EQUAL(v12[jj],v11[jj]);
        }
    }
    
    // Now test everything again using a vector reference
    {
        keycpp::vector_k<int> v1(10,100);
        keycpp::vector_k<int> v2(v1.begin(),v1.size(),1);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        int* it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        BOOST_CHECK_EQUAL(v1[0],1234);
        BOOST_CHECK_EQUAL(v1[9],4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 5678;
        v2[9] = 8765;
        BOOST_CHECK_EQUAL(*v2.begin(),5678);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),8765);
        
        v2.resize(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),20);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        const int* it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        const int* it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
    }
    
    // Now test everything again using a vector reference and an increment
    {
        keycpp::vector_k<int> v1(20,100);
        keycpp::vector_k<int> v2(v1.begin(),10,2);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        int* it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        BOOST_CHECK_EQUAL(v1[0],1234);
        BOOST_CHECK_EQUAL(v1[18],4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 5678;
        v2[9] = 8765;
        BOOST_CHECK_EQUAL(*v2.begin(),5678);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),8765);
        
        v2.resize(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),20);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != nullptr);
        BOOST_CHECK(v2.end() != nullptr);
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        const int* it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        const int* it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
    }
}

// Test KeyCpp Functions
BOOST_AUTO_TEST_CASE(keycpp_test)
{
    // Test fft
    size_t L = 1000;
    double Fs = 1000;
    double T = 1.0/Fs;
    keycpp::vector_k<double> f = Fs*linspace(0.0,1.0,L);
    
    keycpp::vector_k<double> t = keycpp::linspace(0.0,(double)(L-1),L)*T;
    keycpp::vector_k<double> y = sin(2.0*PI*t);
    keycpp::vector_k<std::complex<double>> Y = fft(y)*2.0/1000.0;
    
    BOOST_CHECK(abs(Y[1] - 1.0) < 0.2);
    
}*/
