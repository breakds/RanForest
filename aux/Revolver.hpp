/*********************************************************************************
 * File: Revolver.hpp
 * Description: 1) A random number pool maintainer, maintans groups of random number pool  
 *              2) Capable of group of sample generator without replacement
 * by BreakDS, @ University of Wisconsin-Madison, Apr 6 2012
 * ======================================================================
 *********************************************************************************/


#pragma once

#include <cstdlib>
#include <vector>
#include <limits.h>

using std::vector;

#define SHUFFLER_ERROR UINT_MAX

namespace ran_forest {

  class Shuffler
  {
  private:
    int size;
    vector<int> pool;
    int curPos;
    // pool[address] = id
  public:

    inline Shuffler()
    {
      size = 0;
      pool.clear();
    }
  
    inline Shuffler( const vector<int> &init ) 
    {
      pool.clear();
      size = init.size();
      for ( auto& ele : init ) {
        pool.push_back( ele );
      }
    }

    inline Shuffler( Shuffler &&other ) 
    {
      size = other.size;
      pool.swap( other.pool );
    }

    inline void Reset( const vector<int> &init ) 
    {
      pool.clear();
      size = init.size();
      for ( auto& ele : init ) {
        pool.push_back( ele );
      }
    }

  
    inline void Reset( int n )
    {
      size = n;
      pool.resize( n );
      for ( int i=0; i<n; i++ ) {
        pool[i] = i;
      }
    }
  
    inline Shuffler( int n ) 
    {
      size = n;
      pool.resize( n );
      for ( int i=0; i<n; i++ ) {
        pool[i] = i;
      }
    }

    inline Shuffler( const Shuffler &sfl ) 
    {
      size = sfl.size;
      pool.resize( sfl.size );
      for ( int i=0; i<sfl.size; i++ ) {
        pool[i] = sfl(i);
      }
    }

    inline void operator=( const Shuffler &sfl )
    {
      size = sfl.size;
      pool.resize( sfl.size );
      for ( int i=0; i<sfl.size; i++ ) {
        pool[i] = sfl(i);
      }
    }
  
    inline int Number() const
    {
      return size;
    }
  
    // splice by specifying the id
    // involves looking up
    inline void SpliceID( int id )
    {
      for ( int i=0; i<size-1; i++ ) {
        if ( id == pool[i] ) {
          pool[i] = pool[size-1];
          pool[size-1] = id;
          break;
        }
      }
      size--;
    }
    
    // splice by specifying the address
    // no looking up involved
    inline void SpliceAddress( int add )
    {
      int t = pool[add];
      pool[add] = pool[size-1];
      pool[size-1] = t;
      size--;
    }
  
    inline void Shuffle( int k ) 
    {
      for ( int i=0, end= k < size ? k : size; i < end; i++  ) {
        int r = rand() % (size-i) + i;
        int t = pool[r];
        pool[r] = pool[i];
        pool[i] = t;
      }
    }

    inline void Shuffle() 
    {
      for ( int i=0; i<size; i++  ) {
        int r = rand() % (size-i) + i;
        int t = pool[r];
        pool[r] = pool[i];
        pool[i] = t;
      }
    }

    inline void Keep( int k ) {
      if ( k < size ) {
        size = k;
      }
    }
  
    // Fisher-Yates shuffle
    // - Start Shuffle
    inline void ResetShuffle() {
      curPos = -1;
    }
  
    // - Get Next sample
    inline int Next()
    {
      curPos++;
      if ( curPos >= size ) return SHUFFLER_ERROR; // Out of bound
      int r = rand() % (size-curPos) + curPos;
      int t = pool[r];
      pool[r] = pool[curPos];
      pool[curPos] = t;
      return t;
    }

    // - Delete sample at current position
    inline void Disqualify()
    {
      SpliceAddress( curPos );
      curPos--;
    }




    inline int operator() ( int add ) const
    {
      return pool[add];
    }

    inline void show() const
    {
      printf( "( " );
      for ( int i=0; i<size; i++ ) {
        printf( "%u ", pool[i] );
      }
      printf( ") size = %u\n", size );
    }
  
  
  };

}


