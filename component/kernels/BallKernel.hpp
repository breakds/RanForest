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
#include "LLPack/algorithms/random.hpp"
#include "LLPack/algorithms/algebra.hpp"
#include "LLPack/algorithms/sort.hpp"
#include "LLPack/utils/candy.hpp"


namespace ran_forest
{
  template <typename feature_t, template<typename> class splitter>
  class BallKernel
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
      int depth;

      State( int *i, int l, const Options __attribute__((__unused__)) &options )
        : idx(i), len(l), depth(0) {}

    
      State( int *i, int l, const State& other )
        : idx(i), len(l), depth( other.depth + 1 ) {}

      State( State&& other )
      {
        idx = other.idx;
        len = other.len;
        depth = other.depth;
      }
    };

    const std::vector<feature_t> &dataPoints;
    Options options;

    BallKernel( const std::vector<feature_t>& d, Options o ) : dataPoints(d), options(o) {}
    
    std::vector<int> split( State& state, splitter<dataType> &judger )
    {
      /* exception code:
       * -1 = too few patches within a node
       * -4, -5 = invalid split (totally unbalanced split)
       * -6 = max depth reached
       */

      
      std::vector<int> partition(3);
      partition[0] = 0;
      partition[2] = state.len;
      
      
      if ( state.len <= options.stopNum ) {
        partition[0] = -1;
        return partition;
      }
      
      if ( options.maxDepth == state.depth ) {
        partition[0] = -6;
        return partition;
      }

      
      std::vector<int> hypothesis = rndgen::randperm( state.len, options.numHypo );
      {
        std::vector<double> distances( state.len );
        double bestScore = -1.0;
        double th = 0.0;
        int selected = 0;
        for ( auto& ele : hypothesis ) {
          const feature_t &vantage = dataPoints[state.idx[ele]];
          // calculate distance
          for ( int i=0; i<state.len; i++ ) {
            distances[i] = algebra::dist_l1( vantage, dataPoints[state.idx[i]], options.dim );
          }
          double median = sorting::median( distances );
          // calculate score
          double score = 0.0;
          for ( int i=0; i<state.len; i++ ) {
            score += fabs( distances[i] - median );
          }
          if ( score > bestScore ) {
            bestScore = score;
            th = median;
            selected = state.idx[ele];
          }
        }
        // copy judger
        judger.th = th;
        judger.vantage.resize( options.dim );
        algebra::copy( judger.vantage, dataPoints[selected], options.dim );
      }

      if ( 0 == hypothesis.size() ) {
        partition[0] = -3;
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

  
