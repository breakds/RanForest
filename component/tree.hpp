/* ------------------------------------------------------------------------------------------+
 * This is free and unencumbered software released into the public domain (Unlicense)        |
 *                                                                                           |
 * Source: tree.hpp                                                                          |
 * Author: BreakDS                                                                           |
 * Date: Mon Feb  4 13:19:56 CST 2013                                                        |
 * Description: provide the tree class and its related operations                            |
 * ------------------------------------------------------------------------------------------+
 */

#pragma once
#include <memory>
#include <deque>
#include <string>
#include "LLPack/utils/candy.hpp"
#include "LLPack/utils/extio.hpp"

namespace ran_forest
{

  /* Enum for splitting order */
  enum SplittingOrder {BFS,DFS}; // Breadth First Splitting, Depth First Splitting


  
  
  template <typename dataType, template<typename> class splitter = BinaryOnAxis >
  class Tree
  {
  public:
    // TODO: Simplify kernel options
    class NodeInfo
    {
    public:
      Tree<dataType,splitter> *node;
      int leafID;

      NodeInfo() : node(nullptr), leafID(-1) {}

      inline int GetNodeID()
      {
        return node->nodeID;
      }
    };

    class LeafInfo
    {
    private:
      std::vector<int> store;
    public:
      LeafInfo() : store(0) {}

      inline const int& operator[]( int i )
      {
        return store[i];
      }

      inline void push_back( int e )
      {
        store.push_back( e );
      }
    };

    
  private:
    std::vector<std::unique_ptr< Tree<dataType,splitter> > > child;
    splitter<dataType> judger;
    
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
    template <SplittingOrder order, typename T>
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
    inline void pop( std::deque<T> &q, ENABLE_IF(DFS==order) )
    {
      q.pop_back();
    }
                       

  public:
    /* ---------- constructors ---------- */
    Tree()
    {
      child.clear();
      nodeID = -1;
    }

    template <template <typename,template <typename> class> class kernel, SplittingOrder order = DFS, typename feature_t>
    void grow( const std::vector<feature_t>& dataPoints, // data points
          std::vector<int> &idx,
          std::vector<NodeInfo> &nodes,
          std::vector<LeafInfo> &leaves,
          typename kernel<feature_t,splitter>::Options options )
    {

      typedef kernel<feature_t,splitter> kernelType;
      // TODO: static_assert feature_t's element should have the same type as dataType
      std::deque<std::pair<Tree<dataType,splitter>*,typename kernelType::State> > working_list;

      kernelType core( dataPoints, options );
      
      // push the root into the queue
      working_list.push_back( std::make_pair( this,
                                              typename kernelType::State( &idx[0],
                                                                          static_cast<int>( options.dim ),
                                                                          options.dim ) ) );

      // initialize containers
      nodes.clear();
      leaves.clear();
      
      while ( !working_list.empty() ) {
        Tree<dataType,splitter> *node = fetch<order>(working_list).first;
        typename kernelType::State &state = fetch<order>(working_list).second;
        pop<order>( working_list );


        // emplace the current node
        node->nodeID = static_cast<int>( nodes.size() );
        nodes.emplace( nodes.end() );
        nodes.back().node = this;
        
        // split

        // debugging:
        printf( "---------- nodeID = %d ----------\n", node->nodeID );
        state.shuffler.show();
        std::vector<int> partition = std::move( core.split( state, node->judger ) );
        state.shuffler.show();


        if ( 0 == partition[0] ) {
          /* internal node */
          int branches = static_cast<int>( partition.size() )-1;
          node->child.resize(branches);
          for ( int i=0; i<branches; i++ ) {
            node->child[i].reset( new Tree() );
            working_list.push_back( std::make_pair( node->child[i].get(),
                                                    typename kernelType::State( state.idx + partition[i],
                                                                                partition[i+1] - partition[i],
                                                                                state.shuffler,
                                                                                state.depth + 1 ) ) );
            printf( "len: %d\n", working_list.back().second.len );
            working_list.back().second.shuffler.show();
          }
        } else {
          /* leaf node */
          nodes.back().leafID = leaves.size();
          leaves.emplace( leaves.end() );
          for ( int i=0; i<state.len; i++ ) {
            leaves.back().push_back( idx[i] );
          }
        }

        // debugging:
        char ch;
        scanf( "%c", &ch );

      }
    }

  private:
    /* ---------- private I/Os ---------- */
    void write( FILE* out )
    {
      judger.write( out );
      fwrite( &nodeID, sizeof(int), 1, out );
      unsigned char num = static_cast<unsigned char>( child.size() );
      fwrite( &num, sizeof(unsigned char), 1, out );
      for ( auto& every : child ) {
        every->write( out );
      }
    }

    Tree( FILE *in )
    {
      judger.read( in );
      fread( &nodeID, sizeof(int), 1, in );
      unsigned char num = 0;
      fread( &num, sizeof(unsigned char), 1, in );
      if ( 0 < num ) {
        child.resize( num );
        for ( int i=0; i<num; i++ ) {
          child[i].reset( new Tree( in ) );
        }
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

    static std::unique_ptr<Tree<dataType,splitter> > read( std::string filename )
    {
      std::unique_ptr<Tree<dataType,splitter> > tree;
      WITH_OPEN( in, filename.c_str(), "r" );
      tree.reset( new Tree( in ) );
      END_WITH( in );
      return tree;
    }


    /* ---------- Properties ---------- */
    inline bool isLeaf() const
    {
      return 0 == child.size();
    }

    /* ---------- Queries ---------- */
  };
}
