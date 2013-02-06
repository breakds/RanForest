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

    template <typename T>
    inline int operator()( T p ) const
    {
      // TODO: check element type of p
      if ( p[component] < th ) return 0;
      return 1;
    }
  };
}
