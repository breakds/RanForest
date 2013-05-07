// This file is part of RanForest, a lightweight C++ template library
// for random forest.
//
// By BreakDS <breakds@cs.wisc.edu> - http://www.unlicense.org/ (public domain)
// 
// This file implements the Vantage Point Tree Kernel (packaged as
// class VP).
//
//
// 
// +--- Each kernel class should ---+
//
// 1. have a template paramter for the type of the elements of feature
// vectors that can be handled, called dataType
// ----------------------------------------------------------------------
// 2. specify the splitter that is going to be used and typedef it as
// @typename splitter.
// ----------------------------------------------------------------------
// 3. have a struct named Options. Member variables are options and
// parameters that are used in the ElectSplitter function (described
// below)
// ----------------------------------------------------------------------
// 4. have a class named State. State should contain at least 3 public
// members that describe the state of the current node being splitted:
// - size_t* idx, the array of all the indices currently held at this
//   node
// - size_t len, the length of the array @var idx
// - int depth, the depth of the current node, 0 for root
// ----------------------------------------------------------------------
// 5. have a static function named ElectSplitter. It will generate
// hypothesis for split and pick the best one based on certain
// criteria, as specified by the implementor. This function takes 5
// parameters:
// - a vector of feature_t that contains all the feature vectors used
//   to build the tree.
// - an integer (int) specifies the dimension of each feature vector
// - a @typename State variable ref that describe the current node,
//   see above for details
// - a @typename splitter variable ref that will be used to hold the
//   elected splitter if the election trial is successful
// - a @typename Option variable ref that provides the options
// 
// This function returns a ElectionStatus variable. It will be SUCCESS
// if the election trial is successful, or others as defined in
// tree/define.hpp.



#pragma once

#include <algorithm>
#include "LLPack/algorithms/random.hpp"
#include "../tree/define.hpp"
#include "../splitters/BinaryOnDistance.hpp"

namespace ran_forest
{
  // 1. template paramter 
  template <typename dataType>
  class VP
  {
  public:
    // 2. splitter
    typedef BinaryOnDistance<dataType> splitter;

    // 3. Options
    struct Options
    {
      // stop criteria
      int maxDepth;
      size_t stopNum;
      dataType converge;
      double proportion;
      // other
      size_t numHypo;
      
      Options() : 
        maxDepth(-1), stopNum(5), converge(0), proportion(1.1), numHypo(10) {}
    };

    
    // 4. State
    class State
    {
    public:
      // default properties 
      size_t *idx;
      size_t len;
      int depth;
      
      // extra properties
      
      State( size_t *i, size_t l, const Options __attribute__((__unused__)) &options )
        : idx(i), len(l), depth(0) {}
      
      State( size_t *i, size_t l, const State& other )
        : idx(i), len(l), depth( other.depth + 1 ) {}

      State( State&& other )
      {
        idx = other.idx;
        len = other.len;
        depth = other.depth;
      }
    };


    // 5. Candidate Generator and Elector
    template <typename feature_t>
    static inline ElectionStatus ElectSplitter( const std::vector<feature_t>& dataPoints,
                                                int dim,
                                                State& state,
                                                splitter& judger, 
                                                Options& options )
    {

      if ( state.len < options.stopNum ) {
        return NODE_SIZE_LIMIT_REACHED;
      }

      if ( options.maxDepth == state.depth ) {
        return MAX_DEPTH_REACHED;
      }
      
      std::vector<size_t> vpid = rndgen::randperm( state.len, options.numHypo ); // TODO: size_t of randperm
      double bestScore = -1.0;
      double th = 0.0;
      size_t selected = 0;
      std::vector<double> distances( state.len );
      for ( auto& ele : vpid ) {
        const feature_t &vp = dataPoints[state.idx[ele]];
        for ( size_t i=0; i<state.len; i++ ) {
          distances[i] = algebra::dist_l1( vp, dataPoints[state.idx[i]], dim );
        }
        
        double maxDist = *std::max_element( distances.begin(), distances.end() );
        
        if ( maxDist < options.converge ) {
          return CONVERGED;
        }
        // get median
        std::nth_element( distances.begin(), 
                          distances.begin() + state.len / 2, 
                          distances.end() );
        double median = distances[ state.len / 2 ];
        // calculate score
        for ( size_t i=0; i<state.len; i++ ) {
          distances[i] = fabs( distances[i] - median );
        }
        std::nth_element( distances.begin(),
                          distances.begin() + state.len / 2,
                          distances.end() );
        double score = distances[ state.len / 2 ];

        if ( score > bestScore ) {
          bestScore = score;
          th = median;
          selected = state.idx[ele];
        }
      }
      
      judger.th = th;
      judger.vantage.resize( dim );
      for ( int j=0; j<dim; j++ ) {
        judger.vantage[j] = dataPoints[selected][j];
      }

      return SUCCESS;
    }

  };
  
}
