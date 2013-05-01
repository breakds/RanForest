/*********************************************************************************
 * File: Bipartite.hpp
 * Description: Create a undirected bi-partition graph with weighted edges
 * by BreakDS, University of Wisconsin Madison, Tue Dec 11 16:12:43 CST 2012
 *********************************************************************************/

#pragma once

#include <vector>
#include <utility>
#include <algorithm>

using std::vector;
using std::pair;
using std::make_pair;

namespace ran_forest
{
  class Bipartite
  {
  private:
    vector< vector< pair< size_t, double > > > a_to_b;
    vector< vector< pair< size_t, double > > > b_to_a;

  public:
    Bipartite( size_t numA, size_t numB )
    {
      a_to_b.resize( numA );
      for ( size_t i=0; i<numA; i++ ) {
        a_to_b[i].clear();
      }

      b_to_a.resize( numB );
      for ( size_t i=0; i<numB; i++ ) {
        b_to_a[i].clear();
      }
    }

    Bipartite( Bipartite &&other )
    {
      a_to_b.swap( other.a_to_b );
      b_to_a.swap( other.b_to_a );
    }

    explicit Bipartite( std::string filename )
    {
      WITH_OPEN( in, filename.c_str(), "r" );
      size_t numA = 0;
      size_t numB = 0;
      fread( &numA, sizeof(size_t), 1, in );
      a_to_b.resize( numA );
      fread( &numB, sizeof(size_t), 1, in );
      b_to_a.resize( numB );
      clear();


      for ( size_t a=0; a<numA; a++ ) {
        size_t num = 0;
        fread( &num, sizeof(size_t), 1, in );
        for ( size_t i=0; i<num; i++ ) {
          size_t b = 0;
          double wt = 0.0;
          fread( &b, sizeof(size_t), 1, in );
          fread( &wt, sizeof(double), 1, in );
          add( a, b, wt );
        }
      }
      END_WITH( in );
    }

    void write( std::string filename )
    {
      WITH_OPEN( out, filename.c_str(), "w" );
      size_t numA = sizeA();
      fwrite( &numA, sizeof(int), 1, out );
      size_t numB = sizeB();
      fwrite( &numB, sizeof(int), 1, out );
      for ( size_t a=0; a<numA; a++ ) {
        size_t num = a_to_b[a].size();
        fwrite( &num, sizeof(size_t), 1, out );
        for ( auto& b : a_to_b[a] ) {
          fwrite( &b.first, sizeof(size_t), 1, out );
          fwrite( &b.second, sizeof(double), 1, out );
        }
      }
      END_WITH( out );
    }
    

    const Bipartite& operator=( Bipartite &&other )
    {
      a_to_b.swap( other.a_to_b );
      b_to_a.swap( other.b_to_a );
      return *(this);
    }

    inline void clear()
    {
      for ( auto &ele : a_to_b ) {
        ele.clear();
      }
      for ( auto &ele : b_to_a ) {
        ele.clear();
      }
    }

    inline void grow_a( size_t target ) {
      size_t old = a_to_b.size();
      a_to_b.reserve( target );
      for ( size_t i=old; i<target; i++ ) {
        a_to_b.push_back( vector< pair<size_t, double > >() );
        a_to_b.back().clear();
      }
    }

    inline void grow_b( size_t target ) {
      size_t old = b_to_a.size();
      b_to_a.reserve( target );
      for ( size_t i=old; i<target; i++ ) {
        b_to_a.push_back( vector< pair<size_t, double > >() );
        b_to_a.back().clear();
      }
    }

    
    inline void add( size_t a, size_t b, double wt )
    {
      if ( a >= a_to_b.size() ) {
        grow_a( a + 1 );
      }
      if ( b >= b_to_a.size() ) {
        grow_b( b + 1 );
      }
      a_to_b[a].push_back( make_pair( b, wt ) );
      b_to_a[b].push_back( make_pair( a, wt ) );
    }

    inline vector<pair<size_t,double> >& getSetFrom( size_t a )
    {
      return a_to_b[a];
    }

    inline vector<pair<size_t,double> >& getSetTo( size_t b )
    {
      return b_to_a[b];
    }

    inline const vector<pair<size_t,double> >& from( size_t a ) const
    {
      return a_to_b[a];
    }

    inline const vector<pair<size_t,double> >& to( size_t b ) const
    {
      return b_to_a[b];
    }


    // +-------------------------------------------------------------------------------
    // Property Accessors

    // return the size of set A
    inline size_t sizeA() const 
    {
      return a_to_b.size();
    }

    // return the size of set B
    inline size_t sizeB() const 
    {
      return b_to_a.size();
    }
  };
}
