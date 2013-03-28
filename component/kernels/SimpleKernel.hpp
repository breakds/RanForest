/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: SimpleKernel.hpp                                                                  |
 * Author: BreakDS                                                                           |
 * Date: Tue Feb  5 13:42:53 CST 2013                                                        |
 * Description: This file implement the simple kernel for splitting a tree node              |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include "../Aux/Revolver.hpp"
#include "LLPack/utils/candy.hpp"

namespace ran_forest
{
  template <typename feature_t, template<typename> class splitter>
  class SimpleKernel
  {
  public:
    typedef typename ElementOf<feature_t>::type dataType;

    struct Options
    {
      int maxDepth; // max depth for the forest
      int dim; // dimension of feature
      int numHypo; // number of hypothesis
      int stopNum; // stop splitting when a node contains elements less then stopNum
      dataType converge; // criteria for convergence
      int projDim; // the number of components that will be
                   // considerered in the split process, -1 indicates using all

      Options() : maxDepth(-1), dim(1), numHypo(30), stopNum(5), converge(0), projDim(-1) {}
    }; 
    
    class State
    {
    public:
      int *idx;
      int len;
      Shuffler shuffler;
      int depth;

      State( int *i, int l, int s )
        : idx(i), len(l), shuffler(s), depth(0) {}
      
      State( int *i, int l, const Shuffler& s, int d )
        : idx(i), len(l), shuffler(s), depth(d) {}

      State( State&& other )
      {
        idx = other.idx;
        len = other.len;
        depth = other.depth;
        shuffler = std::move( other.shuffler );
      }
    };

    const std::vector<feature_t> &dataPoints;
    Options options;

    SimpleKernel( const std::vector<feature_t>& d, Options o ) : dataPoints(d), options(o) {}
    
    std::vector<int> split( State& state, splitter<dataType> &judger )
    {
      /* exception code:
       * -1 = too few patches within a node
       * -2 = no candidate component/dimension
       * -3 = all candidate dimension converges
       * -4, -5 = invalid split (totally unbalanced split)
       * -6 = max depth reached
       */

      // if there are more than projDim components available, splice
      // them
      if ( 0 < options.projDim && state.shuffler.Number() > options.projDim ) {
        state.shuffler.ResetShuffle();
        while ( state.shuffler.Number() > options.projDim ) {
          state.shuffler.Next();
          state.shuffler.Disqualify();
        }
      }

      std::vector<int> partition(3);
      partition[0] = 0;
      partition[2] = state.len;
      
      
      if ( state.len <= options.stopNum ) {
        partition[0] = -1;
        return partition;
      }
      if ( 0 == state.shuffler.Number() ) {
        partition[0] = -2;
        return partition;
      }

      if ( options.maxDepth == state.depth ) {
        partition[0] = -6;
        return partition;
      }

      state.shuffler.ResetShuffle();
      
      int trial = 0;
      uint c[options.numHypo];
      dataType th[options.numHypo];
      
      while ( SHUFFLER_ERROR != ( c[trial] = state.shuffler.Next() ) && trial < options.numHypo ) {
        typename Generalized<dataType>::type min = dataPoints[state.idx[0]][c[trial]];
        typename Generalized<dataType>::type max = dataPoints[state.idx[0]][c[trial]];
        for ( int i=1; i<state.len; i++ ) {
          if ( dataPoints[state.idx[i]][c[trial]] > max ) {
            max = dataPoints[state.idx[i]][c[trial]];
          } else if ( dataPoints [state.idx[i]][c[trial]] < min ) {
            min = dataPoints[state.idx[i]][c[trial]];
          }
        }

        if ( max - min < options.converge ) {
          state.shuffler.Disqualify();
        } else {
          dataType range = max - min;
          th[trial] = rand() / static_cast<dataType>( RAND_MAX ) * range * 0.95 + range * 0.025 + min;
          trial++;
        }
      }
      
      if ( 0 == trial ) {
        partition[0] = -3;
        return partition;
      }

      int minDiff = -1;
      for ( int t=0; t<trial; t++ ) {
        int leftNum = 0;
        int rightNum = 0;
        for ( int i=0; i<state.len; i++ ) {
          if ( dataPoints[state.idx[i]][c[t]] < th[t] ) {
            leftNum++;
          } else {
            rightNum++;
          }
        }
        if ( -1 == minDiff || abs( leftNum - rightNum ) < minDiff ) {
          minDiff = abs( leftNum - rightNum );
          judger.th = th[t];
          judger.component = c[t];
        }
      }

      if ( -1 == minDiff ) {
        partition[0] = -4;
        return partition;
      }

      if ( state.len == minDiff ) {
        partition[0] = -5;
        return partition;
      }
      
      int right = -1;
      for ( int i=0; i<state.len; i++ ) {
        if ( 0 == judger( dataPoints[state.idx[i]] ) ) {
          right++;
          int tmp = state.idx[right];
          state.idx[right] = state.idx[i];
          state.idx[i] = tmp;
        }
      }

      partition[1] = right + 1;
      return partition;
    }
  };
}  

  
