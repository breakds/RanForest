/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: forest.hpp                                                                        |
 * Author: BreakDS                                                                           |
 * Date: Wed Feb  6 15:52:28 CST 2013                                                        |
 * Description: a wrapper on trees that maintains a set of trees                             |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include "LLPack/algorithms/random.hpp"
#include "kernels/SimpleKernel.hpp"
#include "splitters/BinaryOnAxis.hpp"
#include "tree.hpp"

namespace ran_forest
{
  template <typename dataType = float, template<typename> class splitter = BinaryOnAxis >
  class Forest
  {

  public:
    typedef Tree<dataType,BinaryOnAxis> treeType;
    typedef typename Tree<dataType,BinaryOnAxis>::NodeInfo NodeInfo;

  private: 
    std::vector<std::unique_ptr<treeType> > trees;
    std::vector<NodeInfo> nodes;
    
  public:

    Forest() : trees(), nodes() {}
    
    // Forest Construction from scratch
    template <template <typename,template <typename> class> class kernel,
              SplittingOrder order = DFS,
              typename feature_t>
    void grow( int n,
               const std::vector<feature_t>& dataPoints,
               typename kernel<feature_t,splitter>::Options options,
               float proportion = 1.1f )

     {
       static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                      "element of feature_t should have the same type as dataType." );
       typedef kernel<feature_t,splitter> kernelType;
       trees.resize( n );

       int len = static_cast<int>( dataPoints.size() );
       int lenPerTree = len;
       if ( proportion < 1.0f ) lenPerTree = static_cast<int>( len * proportion );

       std::vector< std::vector<int> > idx(n);
       std::vector< std::vector<NodeInfo> > tmpNodes(n);

       // construct individual trees
       for ( int i=0; i<n; i++ ) {
         idx[i] = rndgen::randperm( len, lenPerTree );
         trees[i].reset( new treeType );
         trees[i]->grow<kernel,order>( dataPoints,
                                       idx[i],
                                       tmpNodes[i],
                                       options );
       }

       // Merge Node data
       int nodeCount = 0;
       for ( int i=0; i<n; i++ ) {
         nodeCount += static_cast<int>( tmpNodes[i].size() );
       }
       nodes.clear();
       nodes.reserve( nodeCount );
       for ( int i=0; i<n; i++ ) {
         for ( auto& ele : tmpNodes[i] ) {
           ele.node->nodeID = static_cast<int>( nodes.size() );
           nodes.push_back( NodeInfo( std::move( ele ) ) );
         }
       }
     }

     /* ---------- Accessors ---------- */

     inline const NodeInfo& operator[]( int nodeID ) const
     {
       return nodes[nodeID];
     }

     inline int size() const
     {
       return static_cast<int>( trees.size() );
     }

    inline int nodeNum() const
    {
      return static_cast<int>( nodes.size() );
    }

    /* ---------- iterator ---------- */
    typename std::vector<NodeInfo>::const_iterator begin() const
    {
      return nodes.begin();
    }

    typename std::vector<NodeInfo>::const_iterator end() const
    {
      return nodes.end();
    }

    
    /* ---------- Queries ---------- */
    template <typename feature_t>
    inline std::vector<int> query( const feature_t p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      std::vector<int> res;
      res.reserve( trees.size() );

      for ( auto& tree : trees ) {
        res.push_back( tree->query( p ) );
      }

      return res;
    }

  };
}


