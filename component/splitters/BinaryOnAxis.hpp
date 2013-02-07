#pragma once

namespace ran_forest
{
  template <typename dataType>
  class BinaryOnAxis
  {
  public:
    dataType th;
    int component;

    inline void write( FILE *out )
    {
      fwrite( &th, sizeof(dataType), 1, out );
      fwrite( &component, sizeof(int), 1, out );
    }

    inline void read( FILE *in )
    {
      fread( &th, sizeof(dataType), 1, in );
      fread( &component, sizeof(int), 1, in );
    }

    template <typename feature_t>
    inline int operator()( const feature_t& p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      if ( p[component] < th ) return 0;
      return 1;
    }
  };
}
