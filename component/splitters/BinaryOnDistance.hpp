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
  };
  template <typename dataType>
  const std::string BinaryOnDistance<dataType>::name = "Binary On Distance";
}
