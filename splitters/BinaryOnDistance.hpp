// This file is part of RanForest, a lightweight C++ template library
// for random forest.
//
// By BreakDS <breakds@cs.wisc.edu> - http://www.unlicense.org/ (public domain)
// 
// This file implements the class BinaryOnDistance, which is a binary
// tree splitter.
// 
//
// +--- Description on splitter ---+
//
// Each tree node should be associated with a splitter (judge). A
// splitter is a functor (class) that accepts a feature vector, and
// returns an integer (int) indicating which branch (child) the
// feature vector should go to.
//
// Each splitter should have a template parameter dataType which is
// the type of the elements of feature vectors that can be handled,
// and also
//
// 1. a static constant string contains the name of the splitter
// 2. a default constructor that takes no argument
// 3. an operator== for equality comparison
// 4. function write() and read() for serialization
// 5. operator() as it is a functor



#pragma once

#include <string>
#include "LLPack/utils/extio.hpp"
#include "LLPack/algorithms/algebra.hpp"

namespace ran_forest
{
  template <typename dataType>
  class BinaryOnDistance
  {
  public:
    static const std::string name;
  public:
    double th;
    std::vector<dataType> vantage;

    // The default constructor
    BinaryOnDistance() : th(0.0), vantage() {}

    // The move assignment, will be used in tree construction
    const BinaryOnDistance<dataType>& operator==( BinaryOnDistance&& other )
    {
      th = other.th;
      vantage.swap( other );
      return *this;
    }

    inline void write( FILE *out ) const
    {
      fwrite( &th, sizeof(double), 1, out );
      writeVector( out, vantage );
    }

    inline void read( FILE *in )
    {
      fread( &th, sizeof(double), 1, in );
      readVector( in, vantage );
    }
    
    template <typename feature_t>
    inline int operator()( const feature_t& p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      double dist = algebra::dist_l1( p, vantage, static_cast<int>( vantage.size() ) );
      if ( dist < th ) return 0;
      return 1;
    }

    inline bool operator==( const BinaryOnDistance<dataType>& other ) const 
    {
      if ( th != other.th ) return false;
      if ( vantage.size() != other.vantage.size() ) return false;
      for ( int i=0; i<static_cast<int>( vantage.size() ); i++ ) {
        if ( vantage[i] != other.vantage[i] ) return false;
      }
      return true;
    }
  };
  template <typename dataType>
  const std::string BinaryOnDistance<dataType>::name = "Binary On Distance";
}
