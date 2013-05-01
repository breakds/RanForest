// This file is part of RanForest, a lightweight C++ template library
// for random forest.
//
// By BreakDS <breakds@cs.wisc.edu> - http://www.unlicense.org/ (public domain)
// 
// This file implements the random tree/forest, provides the forest
// class and all the relation methods and operations, including
// construction, query, I/O and property accessors. All contents are
// presented under the namespace ran_forest.


#pragma once

#include <cassert>
#include <type_traits>
#include <deque>
#include <functional>
#include "kernels/VP.hpp"


namespace ran_forest
{
  // This is the main class for random forest (or random tree given
  // that the forest contains only one tree.
  //
  // ------------------------------------------------------------
  // The two template arguments are:
  // 
  // - dataType, which can be float, double, unsigned char etc. It
  //   determines the element type of the feature vectors that can be
  //   handled.
  //   
  // - kernel, which is a template class as declared in source codes
  //   under kernel/ subdirectory.
  // ------------------------------------------------------------
  // Data Structures:
  // 1. each tree consists of a set of nodes.
  // 2. node IDs are global, meaning it's unqiue across all the trees.
  // 3. node IDs of children of node i stores in child[i]
  // 4. leaf node i has no children, i.e. child[i].size() == 0
  // 5. judge[i] is a functor that accepts a feature vector, and
  // returns which child of node i the vector should go to next.
  // 6. store[i] contains all the data point IDs (used in tree
  // construction) that falls into (leaf) node i. It is empty if node
  // i is not a leaf.
  // 7. roots[t] stores the node ID of the root of tree t.
  // 8. as an aux data structure, level[i] stores the depth of node i,
  // with root having a depth of 0.
  template <typename dataType = float, template <typename> class kernel = VP>
  class Forest
  {
  private:
    // dimension of the feature vectors
    int dim;
    // for the member variables below, see class description above
    std::vector<size_t> roots; 
    std::vector<std::vector<size_t> > child; 
    std::vector<typename kernel<dataType>::splitter> judge;
    std::vector<int> level;
    std::vector<std::vector<size_t> > store; 
    


  private:
    // The set of functions below operate on deques. Breadth First
    // Search (BFS) pushes to the back of the deque and pops from the
    // front of the deque. Depth First Search (DFS) pushes to the back
    // of the deque as well but also pops from the back. Here fetch
    // and pop are specialized for BFS and DFS, respectively.
    template <SplittingOrder order, typename T>
    inline T& fetch( std::deque<T> &q, ENABLE_IF(BFS==order) )
    {
      return q.front();
    }
    template <SplittingOrder order, typename T>
    inline T& fetch( std::deque<T> &q, ENABLE_IF(DFS==order) )
    {
      return q.back();
    }
    template <SplittingOrder order, typename T>
    inline void pop( std::deque<T> &q, ENABLE_IF(BFS==order) )
    {
      q.pop_front();
    }
    template <SplittingOrder order, typename T>
    inline void pop( std::deque<T> &q, ENABLE_IF(DFS==order) )
    {
      q.pop_back();
    }
    
  public:
    
    // the default constructor
    Forest() : dim(0), roots(), child(), judge(), store() {}

    template <SplittingOrder order = DFS,
              typename feature_t>
    void grow( int n,
               const std::vector<feature_t>& dataPoints,
               int dataDim,
               typename kernel<dataType>::Options options )
    {

      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );

      roots.clear();
      child.clear();
      judge.clear();
      level.clear();
      store.clear();
      

      dim = dataDim;

      size_t len = dataPoints.size();
      size_t lenPerTree = len;
      if ( options.proportion < 1.0 ) 
        lenPerTree = static_cast<size_t>( len * options.proportion );


      roots.resize( n );
      std::vector<std::vector<size_t> > idx(n);
      
      
      ProgressBar progressbar;
      progressbar.reset( n );
      int complete = 0;
#     pragma omp parallel for
      for ( int i=0; i<n; i++ ) {
        idx[i] = rndgen::randperm( len, lenPerTree );
        roots[i] = seed( dataPoints, idx[i], options );
#       pragma omp critical
        {
          progressbar.update( ++complete, "Forest Construction" );
        }
      }
    }
    


    /* grow one tree, return the nodeID of the root */
    template <SplittingOrder order = DFS,
              typename feature_t>
    size_t seed( const std::vector<feature_t>& dataPoints,
                 std::vector<size_t> &idx,
                 typename kernel<dataType>::Options options )
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      std::deque<std::pair<size_t,typename kernel<dataType>::State> > worklist;
      
      size_t root = 0;
#     pragma omp critical
      {
        root = child.size();
        child.emplace_back();
        judge.emplace_back();
        store.emplace_back();
        level.emplace_back( 0 );
      }

      worklist.push_back( std::make_pair( root,
                                          typename kernel<dataType>::State( &idx[0],
                                                                            idx.size(),
                                                                            options ) ) );
      std::vector<int> label( dataPoints.size() );


      while ( !worklist.empty() ) {
        size_t nodeID = fetch<order>(worklist).first;
        typename kernel<dataType>::State state = std::move( fetch<order>(worklist).second );
        pop<order>( worklist );
        
        // try split
        typename kernel<dataType>::splitter newJudge;
        ElectionStatus status = kernel<dataType>::ElectSplitter( dataPoints, dim, state, newJudge, options );
        
        if ( SUCCESS == status ) {
          // assign newJudge
#         pragma omp critical
          {
            judge[nodeID] = std::move( newJudge );
          }
          // calculate branch label
          int maxLabel = -1;
          for ( size_t i=0; i<state.len; i++ ) {
            label[state.idx[i]] = judge[nodeID]( dataPoints[state.idx[i]] );
            if ( label[state.idx[i]] > maxLabel ) {
              maxLabel = label[state.idx[i]];
            }
          }

          // in place counting sort (partition)
          std::vector<size_t> count( maxLabel + 1, 0 );
          for ( size_t i=0; i<state.len; i++ ) count[label[state.idx[i]]]++;
          std::vector<size_t> curpos( maxLabel + 1, 0 );
          for ( int k=1; k<=maxLabel; k++ ) curpos[k] = curpos[k-1] + count[k-1];
          std::vector<size_t> partition( maxLabel + 2, 0 );
          for ( int k=1; k<=maxLabel; k++ ) partition[k] = curpos[k];
          partition[maxLabel+1] = state.len;

          for ( int k=0; k<=maxLabel; k++ ) {
            size_t i = curpos[k];
            while ( i < partition[k+1] ) {
              int k1 = label[state.idx[i]];
              if ( k1 != k ) {
                size_t j = curpos[k1]++;
                size_t tmp = state.idx[i];
                state.idx[i] = state.idx[j];
                state.idx[j] = tmp;
              } else {
                i++;
              }
            }
          }

          // split
#         pragma omp critical
          {
            for ( int k=0; k<=maxLabel; k++ ) {
              size_t id = child.size();
              child.emplace_back();
              judge.emplace_back();
              store.emplace_back();
              level.emplace_back( state.depth + 1 );
              child[nodeID].push_back( id );
              worklist.push_back( std::make_pair( id,
                                                  typename kernel<dataType>::State( state.idx + partition[k],
                                                                                    partition[k+1] - partition[k],
                                                                                    state ) ) );
            }
          }
        } else {
          for ( size_t i=0; i<state.len; i++ ) {
            store[nodeID].push_back( state.idx[i] );
          }
        }
      }
      return root;
    }


    // +-------------------------------------------------------------------------------
    // Input and Output Operations

    void write( std::string dir ) const
    {
      system( strf( "mkdir -p %s", dir.c_str() ).c_str() );
      system( strf( "rm -rf %s/*", dir.c_str() ).c_str() );
      
    }

    void writeTree( std::string dir, int treeID ) const
    {
      WITH_OPEN( in, strf( "%s/tree.%d", dir.c_str(), i ).c_str(), "w" );
      
      END_WITH( in );
    }



    // +-------------------------------------------------------------------------------
    // Query Related Operations
    template <typename feature_t>
    size_t queryTree( const feature_t& p, int treeID ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      size_t i = roots[treeID];
      while ( !child[i].empty() ) {
        i = child[i][judge[i](p)];
      }
      return i;
    }

    template <typename feature_t>
    std::vector<size_t> query( const feature_t& p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      std::vector<size_t> re( roots.size() );
      for ( size_t i=0; i<roots.size(); i++ ) {
        re[i] = queryTree( p, i );
      }
      return re;
    }

    template <typename feature_t>
    std::vector<std::vector<size_t> > batchQuery( const std::vector<feature_t> &dataPoints ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      if ( dataPoints[0].size() != dim ) {
        Error( "RanForest: dimension does not match." );
        exit( -1 );
      }
      std::vector< std::vector<size_t> > re( dataPoints.size() );
      for ( size_t i=0; i<dataPoints.size(); i++ ) {
        re[i].swap( query( dataPoints[i] ) );
      }
      return re;
    }


    /* ---------- Accessors ---------- */
    inline const std::vector<size_t> getStore( size_t nodeID ) const
    {
      return store[nodeID];
    }

    // Return the node IDs of the children of the specified node.
    inline const std::vector<size_t> getChildren( size_t nodeID ) const
    {
      return child[nodeID];
    }
      

    // By default (no parameters provided) it returns the depth of the
    // deepest tree in the forest. Otherwise returns the depth of a
    // particular tree, specified by @param treeID,
    inline int depth( int treeID = -1 ) const
    {
      if ( treeID < 0 ) {
        int re = 0;
        for ( size_t tree=0; tree<roots.size(); tree++ ) {
          int d = depth( tree );
          if ( d > re ) re = d;
        }
        return re;
      }

      // if treeID is not -1, it counts the depth of a particular tree
      assert( treeID < numTrees() );
      return fold_postorder<int>( roots[treeID],
                                  []( const size_t __attribute__((__unused__)) nodeID,
                                      const std::vector<int>& results )
                                  {
                                    return (*std::max_element( results.begin(), results.end() )) + 1;
                                  },
                                  []( const size_t __attribute__((__unused__)) nodeID )
                                  {
                                    return 1;
                                  } );
    }

    inline size_t numLeaves() const
    {
      size_t count = 0; 
      for ( auto& ele : child ) count += ele.empty() ? 1 : 0;
      return count;
    }

    // Return the number of nodes at the depth of @param lv
    inline size_t levelSize( int lv ) const
    {
      size_t count = 0;
      for ( size_t i=0; i<child.size(); i++ ) {
        if ( ( level[i] == lv ) || ( level[i] < lv && child[i].empty() ) ) {
          count++;
        }
      }
      return count;
    }

    // Return all the nodes that are at the depth of @param lv
    inline std::vector<size_t> collectLevel( int lv ) const
    {
      std::vector<size_t> re;
      for ( size_t i=0; i<child.size(); i++ ) {
        if ( ( level[i] == lv ) || ( level[i] < lv && child[i].empty() ) ) {
          re.push_back( i );
        }
      }
      return re;
    }

    inline size_t numNodes() const
    {
      return child.size();
    }

    inline int numTrees() const
    {
      return static_cast<int>( roots.size() );
    }

    inline size_t treeRoot( int treeID )
    {
      return roots[treeID];
    }

    
    /* ---------- Extra ---------- */

    // This function applies @param fun at each internal node, and
    // @param base at each leaf node, with a postorder traversal. For
    // each internal node, the return value of @param internal_fun
    // will be based on both the node itself and the return values of
    // either @param internal_fun or @param leaf_fun on its children,
    // as stored in a vector<T>.
    //
    // The first parameter for both internal_fun and leaf_fun is the
    // nodeID.
    // 
    // TODO (BreakDS): This is a recursive implementation. For
    // efficiency a iterative implementation may be better. However
    // this is not of the highest priority since @fun fold_postorder
    // is mainly used for debugging purpose.
    template <typename T>
    T fold_postorder( size_t nodeID,
                      std::function<T(const size_t, const std::vector<T>&)> internal_fun,
                      std::function<T(const size_t)> leaf_fun ) const
    {
      if ( child[nodeID].empty() ) {
        return leaf_fun( nodeID );
      } else {
        std::vector<T> results( child[nodeID].size() );
        for ( size_t c=0; c<child[nodeID].size(); c++ ) {
          results[c] = fold_postorder( child[nodeID][c], internal_fun, leaf_fun );
        }
        return internal_fun( nodeID, results );
      }
    }
  };
  
}
