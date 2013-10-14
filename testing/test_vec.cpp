#include <boost/test/unit_test.hpp>
#include "keycpp.h"

BOOST_AUTO_TEST_SUITE(KeyCpp_Unit_Testing)

BOOST_AUTO_TEST_CASE(vector_k_test)
{
    {
        keycpp::vector_k<int> v1;
        BOOST_CHECK_EQUAL(v1.capacity(),0);
        BOOST_CHECK_EQUAL(v1.size(),0);
        BOOST_CHECK(v1.empty());
        BOOST_CHECK(v1.begin() == v1.end());
        
        keycpp::vector_k<int> v2(10);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        keycpp::vector_k<int>::iterator it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
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
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        keycpp::vector_k<int>::const_iterator it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        keycpp::vector_k<int>::const_iterator it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
        
        v2.push_back(3333);
        BOOST_CHECK(v2.capacity() >= 21);
        BOOST_CHECK_EQUAL(v2.size(),21);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
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
        BOOST_CHECK(v2.begin() != v2.end());
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
        keycpp::vector_k<int> v2(&v1[0],v1.size(),1);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        keycpp::vector_k<int>::iterator it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        BOOST_CHECK_EQUAL(v1[0],1234);
        BOOST_CHECK_EQUAL(v1[9],4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
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
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        keycpp::vector_k<int>::const_iterator it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        keycpp::vector_k<int>::const_iterator it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
    }
    
    // Now test everything again using a vector reference and an increment
    {
        keycpp::vector_k<int> v1(20,100);
        keycpp::vector_k<int> v2(&v1[0],10,2);
        BOOST_CHECK(v2.capacity() >= 10);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1234;
        v2[9] = 4321;
        BOOST_CHECK_EQUAL(*v2.begin(),1234);
        keycpp::vector_k<int>::iterator it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),4321);
        BOOST_CHECK_EQUAL(v1[0],1234);
        BOOST_CHECK_EQUAL(v1[18],4321);
        
        v2.reserve(20);
        BOOST_CHECK(v2.capacity() >= 20);
        BOOST_CHECK_EQUAL(v2.size(),10);
        BOOST_CHECK(!v2.empty());
        BOOST_CHECK(v2.begin() != v2.end());
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
        BOOST_CHECK(v2.begin() != v2.end());
        v2[0] = 1111;
        v2[19] = 2222;
        BOOST_CHECK_EQUAL(*v2.begin(),1111);
        it_end = v2.end();
        std::advance(it_end,-1);
        BOOST_CHECK_EQUAL(*(it_end),2222);
        keycpp::vector_k<int>::const_iterator it_begin2 = v2.begin();
        BOOST_CHECK_EQUAL(*it_begin2,1111);
        keycpp::vector_k<int>::const_iterator it_end2 = v2.end();
        std::advance(it_end2,-1);
        BOOST_CHECK_EQUAL(*(it_end2),2222);
        BOOST_CHECK_EQUAL(v2.front(),1111);
        BOOST_CHECK_EQUAL(v2.back(),2222);
    }
    
    // Test initialization values
    {
        keycpp::vector_k<int> v1(10);
        BOOST_CHECK_EQUAL(v1[0],0);
        BOOST_CHECK_EQUAL(v1[9],0);
        
        double tol = 1e-10;
        keycpp::vector_k<double> v2(10);
        BOOST_CHECK_CLOSE(v2[0],0.0,tol);
        BOOST_CHECK_CLOSE(v2[9],0.0,tol);
    }
    
    // Test vector_k iterator
    {
        keycpp::vector_k<int> v1(10);
        for(keycpp::vector_k<int>::iterator it = v1.begin(); it != v1.end(); it++)
        {
            *it = 1;
        }
        for(auto ii : v1)
        {
            BOOST_CHECK_EQUAL(ii,1);
        }
        for(keycpp::vector_k<int>::iterator it = v1.begin(); it != v1.end(); it += 1)
        {
            BOOST_CHECK_EQUAL(*it,1);
        }
        keycpp::vector_k<int>::iterator it1 = v1.begin() + 5;
        BOOST_CHECK_EQUAL(*it1,1);
        
        keycpp::vector_k<int> v2(&v1[0],5,2);
        for(keycpp::vector_k<int>::iterator it = v2.begin(); it != v2.end(); it++)
        {
            *it = 2;
        }
        for(keycpp::vector_k<int>::iterator it = v2.begin(); it != v2.end(); it++)
        {
            BOOST_CHECK_EQUAL(*it,2);
        }
        
        keycpp::vector_k<int>::iterator it2 = v1.begin() + 2;
        BOOST_CHECK_EQUAL(*it2,2);
        BOOST_CHECK_EQUAL(v1[0],2);
        BOOST_CHECK_EQUAL(v1[1],1);
        
        keycpp::vector_k<int>::iterator it3 = v2.begin() + 1;
        BOOST_CHECK_EQUAL(it3.get_inc(),2);
        BOOST_CHECK_EQUAL(*it3,2);
        it3 += 1;
        BOOST_CHECK_EQUAL(*it3,2);
        
        for(auto ii : v2)
        {
            BOOST_CHECK_EQUAL(ii,2);
        }
        
        keycpp::vector_k<int>::iterator it4 = v2.end();
        BOOST_CHECK_EQUAL((--it4) - v2.begin(),v2.size()-1);
        
        keycpp::vector_k<int>::iterator it5 = v2.begin();
        BOOST_CHECK_EQUAL(v2.end() - (++it5),v2.size()-1);
    }
    
    // Test with Boost.range
    {
        using namespace boost;
        typedef keycpp::vector_k<int> vec_t;
        vec_t                    vec;
        vec.push_back( 3 ); vec.push_back( 4 );
        const vec_t              cvec( vec ); 

        BOOST_STATIC_ASSERT(( is_same< range_value<vec_t>::type, vec_t::value_type >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_iterator<vec_t>::type, vec_t::iterator >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_iterator<const vec_t>::type, vec_t::const_iterator >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_difference<vec_t>::type, vec_t::difference_type >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_size<vec_t>::type, vec_t::size_type >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_iterator<vec_t>::type, vec_t::iterator >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_iterator<const vec_t>::type, vec_t::const_iterator >::value ));
        
        BOOST_STATIC_ASSERT(( is_same< range_difference<const vec_t>::type, vec_t::difference_type >::value ));
        BOOST_STATIC_ASSERT(( is_same< range_size<const vec_t>::type, vec_t::size_type >::value ));

        BOOST_CHECK( begin( vec ) == vec.begin() );
        BOOST_CHECK( end( vec )   == vec.end() );
        BOOST_CHECK( empty( vec ) == vec.empty() );
        BOOST_CHECK( (std::size_t)size( vec ) == vec.size() );
        
        BOOST_CHECK( begin( cvec ) == cvec.begin() );
        BOOST_CHECK( end( cvec )   == cvec.end() );
        BOOST_CHECK( empty( cvec ) == cvec.empty() );
        BOOST_CHECK( (std::size_t)size( cvec ) == cvec.size() );
    }
}

BOOST_AUTO_TEST_SUITE_END()
