#pragma once
#include "define.hpp"

namespace ran_forest
{
  template <typename dataType>
  class VP
  {
  public:
    typedef BinaryOnDistance<dataType> splitter;

    static const int branches;


    // 1. Options 
    struct Options
    {
      // stop criteria
      int maxDepth;
      size_t stopNum;
      dataType converge;
      // other
      size_t numHypo;
      
      Options() : maxDepth(-1), numHypo(10), stopNum(5), converge(0) {}
    };

    
    // 2. State
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
        l = other.l;
        depth = other.depth;
      }
    };


    // 3. Generate Candidate
    template <typename feature_t>
    static inline ElectionStatus ElectSplitter( const std::vector<feature_t>& dataPoints,
                                                int dim,
                                                State& state,
                                                splitter& judger )
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
      for ( auto& ele : vpid ) {
        const feature_t &vp = dataPoints[state.idx[ele]];
        for ( size_t i=0; i<state.len; i++ ) {
          distance[i] = algebra::dist_l1( vp, dataPoints[state.idx[i]], dim );
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
      
      return SUCCESS;
    }

  };
  
}
