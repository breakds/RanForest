/* This is free and unencumbered software released into the public domain (Unlicense)
 *
 * Source: tree.hpp
 * Author: BreakDS
 * Date: Mon Feb  4 13:19:56 CST 2013
 * Description: provide the tree class and its related operations
 */

#pragma once
#include "LLPack/utils/candy.hpp"

namespace ran_forest
{

  /* Enum for splitting order */
  enum SplittingOrder {BFS,DFS}; // Breadth First Splitting, Depth First Splitting
  
  template <typename dataType, template<typename> class splitter = SimpleSplitter >
  class Tree
  {
  private:
    std::unique_ptr< Tree<dataType,splitter> > child;
    splitter judger;
    
  public:
    int nodeID;


  private:
    // fetch next from working list, breadth first version
    template <SplittingOrder order, typename T>
    inline T& fetch( std::deque<T> &q, ENABLE_IF(BFS==order) )
    {
      return q.front();
    }

    // fetch next from working list, depth first version
    inline T& fetch( std::deque<T> &q, ENABLE_IF(DFS==order) )
    {
      return q.back();
    }

    // pop one element from working list, breadth first version
    template <SplittingOrder order, typename T>
    inline void pop( std::deque<T> &q, ENABLE_IF(BFS==order) )
    {
      q.pop_front();
    }

    // pop one element from working list, depth first version
    template <SplittingOrder order, typename T>
    inline void pop( std::deque<T> &q, ENABLE_IF(BFS==order) )
    {
      q.pop_back();
    }
                       

    /* ---------- constructors ---------- */
    Tree()
    {
      child.clear();
      nodeID = -1;
    }

    template <template <typename,typename> class kernel, SplittingOrder order = DFS, typename feature_t>
    Tree ( const std::vector<feature_t>& dataPoints, // data points
           const std::vector<int> &idx,
           std::vector<LeafInfo> &leaves,
           kernel<feature_splitter>::Options options )
    {

      // TODO: static_assert feature_t's element should have the same type as dataType
      std::deque<std::pair<Tree<kernel>*,typename kernel::State> > working_list;

      typename kernel<feature_t,splitter> core( dataPoints, options );

      // push the root into the queue
      working_list.push_back( std::make_pair( this,
                                              typename kernel::State( &idx[0],
                                                                      static_cast<int>( idx.size() ),
                                                                      options.dim ) ) );

      int nodeCount = 0;
      while ( !stack.empty() ) {
        Tree<kernel> *node = fetch<order>(working_list).first;
        typename kernel::State &state = fetch<order>(working_list).second;
        pop<order>( working_list );
        node->nodeID = nodeCount++;

        // split 
        std::vector<int> partition = std::move( core.split( state, node->judger ) );

        if ( 0 == partition[0] ) {
          /* internal node */
          int branches = statitc_cast<int>( partition.size() )-1;
          for ( int i=0; i<branches; i++ ) {
            node->child[i].reset( new Tree() );
            working_list.push_back( std::make_pair( node->child[i].get(),
                                                    typename kernel::State( state.idx + partition[i],
                                                                            partition[i+1] - partition[i],
                                                                            state.shuffler,
                                                                            state.depth + 1 ) ) );
          }
        } else {
          /* leaf node */
          // TODO: push idx information into leaves
        }
      }
    }

  private:
    /* ---------- private I/Os ---------- */
    void write( FILE* out )
    {
      judger.write( out );
      unsigned char num = static_cast<unsigned char>( child.size() );
      fwrite( &num, sizeof(unsigned char), 1, out );
      for ( auto& every : child ) {
        every->write( out );
      }
    }
    
  public:
    /* ---------- public I/Os ---------- */
    void write( std::string filename )
    {
      WITH_OPEN( out, filename.c_str(), "w" );
      write( out );
      END_WITH( out );
    }
  };
}
