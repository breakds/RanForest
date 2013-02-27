/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: MaxGapSubspaceKernel.hpp                                                          |
 * Author: BreakDS                                                                           |
 * Date: Tue Feb  5 13:42:53 CST 2013                                                        |
 * Description: This file implement the kernel for splitting a tree                          |
 * node while the criteria is picking the largest gap on the given axis                      |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include "../Aux/Revolver.hpp"
#include "LLPack/utils/candy.hpp"

namespace ran_forest
{
  template <typename feature_t, template<typename> class splitter>
  class MaxGapSubspaceKernel
  {
  public:
    typedef typename ElementOf<feature_t>::type dataType;

    struct Options
    {
      int maxDepth; // max depth for the forest, -1 for no limit
      int dim; // # dimension of feature
      int dimPrelim; // # dimension of the preliminary axis
      int dimFinal; // # dimension of the final axis
      int numHypo; // # of hypothesis prelimnary axis
      int numHypoTh; // # of hypothesis thresholds for every final axis 
      int stopNum; // stop splitting when a node contains elements less then stopNum
      dataType converge; // criteria for component convergence

      Options() : maxDepth(-1), dim(1), dimPrelim(5), 
                  dimFinal(3), numHypo(1), numHypoTH(10),
                  stopNum(5), converge(0) {}
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

    MaxGapSubspaceKernel( const std::vector<feature_t>& d, Options o ) : dataPoints(d), options(o) {}
    
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


      
      int trial = 0;
      
      // Sample Preliminary Subspace
      state.shuffler.ResetShuffle();
      int realDimPrelim = 0;
      int c[options.dimPrelim];
      while ( SHUFFLER_ERROR != ( c[realDimPrelim] = state.shuffler.Next() ) && realDimPrelim < options.numHypo ) {
        typename Generalized<dataType>::type min = dataPoints[state.idx[0]][c[realDimPrelim]];
        typename Generalized<dataType>::type max = dataPoints[state.idx[0]][c[realDimPrelim]];
        for ( int i=1; i<state.len; i++ ) {
          if ( dataPoints[state.idx[i]][c[realDimPrelim]] > max ) {
            max = dataPoints[state.idx[i]][c[realDimPrelim]];
          } else if ( dataPoints [state.idx[i]][c[realDimPrelim]] < min ) {
            min = dataPoints[state.idx[i]][c[realDimPrelim]];
          }
        }

        if ( max - min < options.converge ) {
          state.shuffler.Disqualify();
        } else {
          realDimPrelim++;
        }
      }
      
      if ( 0 == realDimPrelim ) {
        partition[0] = -3;
        return partition;
      }


      // Generate Random Vector for Projection
      // TODO -> axis[];
      
      // farthest away minIdx, maxIdx

      // 
      

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

  
