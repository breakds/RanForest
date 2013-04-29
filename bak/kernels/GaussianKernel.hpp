/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: GaussianKernel.hpp                                                                |
 * Author: BreakDS                                                                           |
 * Date: Tue Feb  5 13:42:53 CST 2013                                                        |
 * Description: This file implement the distribution criteria for splitting a node           |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include "../Aux/Revolver.hpp"
#include "LLPack/utils/candy.hpp"

namespace ran_forest
{
  template <typename feature_t, template<typename> class splitter>
  class GaussianKernel
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

      Options() : maxDepth(-1), dim(1), numHypo(30), stopNum(5), converge(0) {}
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

    GaussianKernel( const std::vector<feature_t>& d, Options o ) : dataPoints(d), options(o) {}
    
    std::vector<int> split( State& state, splitter<dataType> &judger )
    {
      /* exception code:
       * -1 = too few patches within a node
       * -2 = no candidate component/dimension
       * -3 = all candidate dimension converges
       * -4, -5 = invalid split (totally unbalanced split)
       * -6 = max depth reached
       */

      std::vector<int> partition(3);
      partition[0] = 0;
      partition[2] = state.len;
      
      // too few data points in a node
      if ( state.len <= options.stopNum ) {
        partition[0] = -1;
        return partition;
      }

      // no component ready to be hypothesis
      if ( 0 == state.shuffler.Number() ) {
        partition[0] = -2;
        return partition;
      }

      // reaching max depth limit
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

      double minVariance = -1.0;
      for ( int t=0; t<trial; t++ ) {
        int leftNum = 0;
        int rightNum = 0;
        typename Generalized<dataType>::type leftSum = 0;
        typename Generalized<dataType>::type rightSum = 0;
        typename Generalized<dataType>::type leftMomentum = 0;
        typename Generalized<dataType>::type rightMomentum = 0;
        for ( int i=0; i<state.len; i++ ) {
          const dataType &x = dataPoints[state.idx[i]][c[t]];
          if ( x < th[t] ) {
            leftNum++;
            leftSum += x;
            leftMomentum += x * x;
          } else {
            rightNum++;
            rightSum += x;
            rightMomentum += x * x;
          }
        }
        
        double leftMean(0.0), leftVariance(0.0), rightMean(0.0), rightVariance(0.0);
        if ( leftNum > 0 ) {
          double inv = 1.0 / leftNum;
          leftMean = leftSum * inv;
          leftVariance = leftMomentum * inv - leftMean * leftMean;
        }
        if ( rightNum > 0 ) {
          double inv = 1.0 / rightNum;
          rightMean = rightSum * inv;
          rightVariance = rightMomentum * inv - rightMean * rightMean;
        }

        double variance = leftVariance + rightVariance;
        
        if ( 0 > minVariance || variance < minVariance ) {
          minVariance = variance;
          judger.th = th[t];
          judger.component = c[t];
        }
      }

      if ( 0 > minVariance ) {
        partition[0] = -4;
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

  
