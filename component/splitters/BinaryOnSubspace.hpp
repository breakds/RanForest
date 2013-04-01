#pragma once

#include <string>
#include "LLPack/utils/extio.hpp"
#include "LLPack/utils/candy.hpp"

namespace ran_forest
{
  template <typename dataType>
  class BinaryOnSubspace
  {
  public:
    static const std::string name;
  public:
    dataType th;
    std::vector<int> components;
    std::vector<dataType> projaxis;
    
    inline void write( FILE *out ) const
    {
      fwrite( &th, sizeof(dataType), 1, out );
      writeVector( out, components );
      writeVector( out, projaxis );
    }

    inline void read( FILE *in )
    {
      fread( &th, sizeof(dataType), 1, in );
      readVector( in, components );
      readVector( in, projaxis );
    }

    template <typename feature_t>
    inline int operator()( const feature_t& p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      typename Generalized<dataType>::type val = 0;
      int i = 0;
      for ( auto& ele : components ) {
        val += projaxis[i++] * p[ele];
      }
      if ( val < th ) return 0;
      return 1;
    }
  };
  template <typename dataType>
  const std::string BinaryOnSubspace<dataType>::name = "Binary On Subspace";
}
