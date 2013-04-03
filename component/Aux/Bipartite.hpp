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
    vector< vector< pair< int, double > > > a_to_b;
    vector< vector< pair< int, double > > > b_to_a;

  public:
    Bipartite( int numA, int numB )
    {
      a_to_b.resize( numA );
      for ( int i=0; i<numA; i++ ) {
        a_to_b[i].clear();
      }

      b_to_a.resize( numB );
      for ( int i=0; i<numB; i++ ) {
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
      int numA = 0;
      int numB = 0;
      fread( &numA, sizeof(int), 1, in );
      a_to_b.resize( numA );
      fread( &numB, sizeof(int), 1, in );
      b_to_a.resize( numB );
      clear();


      for ( int a=0; a<numA; a++ ) {
        int num = 0;
        fread( &num, sizeof(int), 1, in );
        for ( int i=0; i<num; i++ ) {
          int b = 0;
          double wt = 0.0;
          fread( &b, sizeof(int), 1, in );
          fread( &wt, sizeof(double), 1, in );
          add( a, b, wt );
        }
      }
      END_WITH( in );
    }

    void write( std::string filename )
    {
      WITH_OPEN( out, filename.c_str(), "w" );
      int numA = sizeA();
      fwrite( &numA, sizeof(int), 1, out );
      int numB = sizeB();
      fwrite( &numB, sizeof(int), 1, out );
      for ( int a=0; a<numA; a++ ) {
        int num = static_cast<int>( a_to_b[a].size() );
        fwrite( &num, sizeof(int), 1, out );
        for ( auto& b : a_to_b[a] ) {
          fwrite( &b.first, sizeof(int), 1, out );
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

    inline void grow_a( int target ) {
      int old = static_cast<int>( a_to_b.size() );
      a_to_b.reserve( target );
      for ( int i=old; i<target; i++ ) {
        a_to_b.push_back( vector< pair<int, double > >() );
        a_to_b.back().clear();
      }
    }

    inline void grow_b( int target ) {
      int old = static_cast<int>( b_to_a.size() );
      b_to_a.reserve( target );
      for ( int i=old; i<target; i++ ) {
        b_to_a.push_back( vector< pair<int, double > >() );
        b_to_a.back().clear();
      }
    }

    
    inline void add( int a, int b, double wt )
    {
      if ( a >= static_cast<int>( a_to_b.size() ) ) {
        grow_a( a + 1 );
      }
      if ( b >= static_cast<int>( b_to_a.size() ) ) {
        grow_b( b + 1 );
      }
      a_to_b[a].push_back( make_pair( b, wt ) );
      b_to_a[b].push_back( make_pair( a, wt ) );
    }

    inline vector<pair<int,double> >& getSetFrom( int a )
    {
      return a_to_b[a];
    }

    inline vector<pair<int,double> >& getSetTo( int b )
    {
      return b_to_a[b];
    }

    inline const vector<pair<int,double> >& from( int a ) const
    {
      return a_to_b[a];
    }

    inline const vector<pair<int,double> >& to( int b ) const
    {
      return b_to_a[b];
    }

    
    inline int sizeA() const 
    {
      return static_cast<int>( a_to_b.size() );
    }

    inline int sizeB() const 
    {
      return static_cast<int>( b_to_a.size() );
    }

    /* ---------- Developer's Operations ---------- */

    size_t MedianContainmentA( const std::vector<int> *idx = nullptr )
    {
      std::vector<int> indices;
      if ( nullptr == idx ) {
        indices.resize( a_to_b.size() );
        for ( size_t i=0; i<a_to_b.size(); i++ ) {
          indices[i] = i;
        }
      } else {
        indices = *idx;
      }
      
      for ( size_t i=0; i<a_to_b.size(); i++ ) {
        indices[i] = i;
      }
      std::nth_element( indices.begin(),
                        indices.begin() + indices.size() / 2,
                        indices.end(),
                        [this]( int i, int j )
                        {
                          return a_to_b[i].size() < a_to_b[j].size();
                        } );
      return a_to_b[indices[indices.size()/2]].size();
    }
    
    size_t MedianContainmentB( const std::vector<int> *idx = nullptr )
    {
      std::vector<int> indices;
      if ( nullptr == idx ) {
        indices.resize( b_to_a.size() );
        for ( size_t i=0; i<b_to_a.size(); i++ ) {
          indices[i] = i;
        }
      } else {
        indices = *idx;
      }
      std::nth_element( indices.begin(),
                        indices.begin() + indices.size() / 2,
                        indices.end(),
                        [this]( int i, int j )
                        {
                          return b_to_a[i].size() < b_to_a[j].size();
                        } );
      return b_to_a[indices[indices.size()/2]].size();
    }
  };
}
