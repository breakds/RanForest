/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: MaxGapSubspaceKernel.hpp                                                          |
 * Author: BreakDS                                                                           |
 * Date: Wed Feb 27 15:26:31 CST 2013                                                        |
 * Description: This file implement the kernel for splitting a tree                          |
 * node while the criteria is picking the largest gap on the given axis                      |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include <random>
#include <chrono>
#include <algorithm>
#include "../Aux/Revolver.hpp"
#include "LLPack/utils/candy.hpp"
#include "LLPack/algorithms/random.hpp"
#include "LLPack/algorithms/sort.hpp"
#include "LLPack/algorithms/algebra.hpp"


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
      int numHypoTH; // # of hypothesis thresholds for every final axis 
      int stopNum; // stop splitting when a node contains elements less then stopNum
      dataType converge; // criteria for component convergence
      float margin;
      std::mt19937 rng;

      Options() : maxDepth(-1), dim(1), dimPrelim(5), 
                  dimFinal(3), numHypo(1), numHypoTH(10),
                  stopNum(5), converge(0), margin(0.1) {
        rng.seed( std::chrono::system_clock::now().time_since_epoch().count() );
      }
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


  private:
    inline double dotprod( const feature_t& feat, const std::vector<int>& comp,
                    const std::vector<double>&axis )
    {
      int i = 0;
      double s = 0.0;
      for ( auto& ele : axis ) {
        s += feat[comp[i++]] * ele;
      }
      return s;
    }
    
  public:
    std::vector<int> split( State& state, splitter<dataType> &judger )
    {
      /* exception code:
       * -1 = too few patches within a node
       * -2 = no candidate component/dimension
       * -3 = all candidate dimension converges
       * -4 = invalid split (totally unbalanced split)
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


      
      double bestScore = -state.len-1;
      for ( int trial = 0; trial < options.numHypo; trial++ ) {
        // Sample Preliminary Subspace
        state.shuffler.ResetShuffle();
        int realDimPrelim = 0;
        std::vector<uint> c0(options.dimPrelim);
      
        while ( realDimPrelim < options.dimPrelim
                && SHUFFLER_ERROR != ( c0[realDimPrelim] = state.shuffler.Next() ) ) {
          typename Generalized<dataType>::type min = dataPoints[state.idx[0]][c0[realDimPrelim]];
          typename Generalized<dataType>::type max = dataPoints[state.idx[0]][c0[realDimPrelim]];
          for ( int i=1; i<state.len; i++ ) {
            if ( dataPoints[state.idx[i]][c0[realDimPrelim]] > max ) {
              max = dataPoints[state.idx[i]][c0[realDimPrelim]];
            } else if ( dataPoints [state.idx[i]][c0[realDimPrelim]] < min ) {
              min = dataPoints[state.idx[i]][c0[realDimPrelim]];
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

        std::vector<int> c( realDimPrelim );
        for ( int i=0; i<realDimPrelim; i++ ) c[i] = static_cast<int>( c0[i] );


      
        // Generate Random Vector for Projection, and find the max/min
        int maxIdx = 0;
        int minIdx = 0;
        double maxPrj = 0.0;
        double minPrj = 0.0;
        std::vector<double> axis;
        {
          axis = rndgen::rnd_unit_vec<double>( realDimPrelim, options.rng );
          maxPrj = dotprod( dataPoints[state.idx[0]], c, axis );
          minPrj = maxPrj;

          for ( int i=1; i<state.len; i++ ) {
            double prj = dotprod( dataPoints[state.idx[i]], c, axis );
            if ( prj > maxPrj ) {
              maxPrj = prj;
              maxIdx = i;
            } else if ( prj < minPrj ) {
              minPrj = prj;
              minIdx = i;
            }
          }
        }

        // find the most significant components between difference of max and min
        // and generate final axis
        c.resize( options.dimFinal );
        axis.resize( options.dimFinal );
        {
          std::vector<double> diff(options.dim);
          for ( int j=0; j<options.dim; j++ ) {
            diff[j] = fabs( dataPoints[state.idx[maxIdx]][j] -
                            dataPoints[state.idx[minIdx]][j] );
          }
          std::vector<int> sorted = sorting::index_sort( diff );

          for ( int j=0; j<options.dimFinal; j++ ) {
            c[j] = sorted[j];
            axis[j] = dataPoints[state.idx[maxIdx]][j] - dataPoints[state.idx[minIdx]][j];
          }
          double len = norm_l2( &axis[0], options.dimFinal );
          for ( auto& ele : axis ) ele /= len;
        }


        // Find min proj and max proj
        maxPrj = dotprod( dataPoints[state.idx[0]], c, axis );
        minPrj = maxPrj;
        for ( int i=1; i<state.len; i++ ) {
          double prj = dotprod( dataPoints[state.idx[i]], c, axis );
          if ( prj > maxPrj ) {
            maxPrj = prj;
          } else if ( prj < minPrj ) {
            minPrj = prj;
          }
        }
      

        // generate candidate for thresholds
        {
          std::vector<dataType> thCands = rndgen::rnd_uniform_real<float>( options.numHypoTH,
                                                                           maxPrj * options.margin +
                                                                           minPrj * ( 1.0f - options.margin ),
                                                                           maxPrj * ( 1.0f - options.margin ) +
                                                                           minPrj * options.margin,
                                                                           options.rng );
          std::vector<int> leftNum( options.numHypoTH );
          std::vector<int> rightNum( options.numHypoTH );
          std::vector<dataType> lowerBound( options.numHypoTH, minPrj );
          std::vector<dataType> upperBound( options.numHypoTH, maxPrj );


          for ( int i=0; i<state.len; i++ ) {
            double prj = dotprod( dataPoints[state.idx[i]], c, axis );
            for ( int j=0; j<options.numHypoTH; j++ ) {
              if ( prj < thCands[j] ) {
                leftNum[j]++;
                if ( prj > lowerBound[j] ) lowerBound[j] = prj;
              } else {
                rightNum[j]++;
                if ( prj < upperBound[j] ) upperBound[j] = prj;
              }
            }
          }

          // find max gap
          int picked = -1;
          double maxGap = -1.0;
          for ( int j=0; j<options.numHypoTH; j++ ) {
            if ( 0 < leftNum[j] && 0 < rightNum[j] ) {
              double gap = upperBound[j] - lowerBound[j];

              if ( -1 == picked || gap > maxGap ) {
                maxGap = gap;
                picked = j;
              }
            }
          }
    
          if ( -1 == picked ) {
            partition[0] = -4;
          }

          // double score = maxGap / ( maxPrj - minPrj );
          double score = -abs( leftNum[picked] - rightNum[picked] );
          // update judger
          if ( score > bestScore ) {
            bestScore = score;
            judger.th = thCands[picked];
            judger.components.swap( c );
            {
              judger.projaxis.resize( axis.size() );
              int j = 0;
              for ( auto& ele : axis ) {
                judger.projaxis[j++] = ele;
              }
            }
          }
        }
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

  
