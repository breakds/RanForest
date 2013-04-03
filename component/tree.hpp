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
      std::vector<int> store;
      
      NodeInfo() : node(nullptr), store() {}

      explicit NodeInfo( FILE *in )
      {
        node = nullptr;
        readVector( in, store );
      }

      explicit NodeInfo( NodeInfo&& other )
      {
        node = other.node;
        store.swap( other.store );
      }

      inline int GetNodeID()
      {
        return node->nodeID;
      }

      inline bool isLeaf() const
      {
        return node->isLeaf();
      }

      inline void write( FILE *out ) const
      {
        writeVector( out, store );
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
               typename kernel<feature_t,splitter>::Options options )
    {


      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      
      typedef kernel<feature_t,splitter> kernelType;
      std::deque<std::pair<Tree<dataType,splitter>*,typename kernelType::State> > working_list;

      kernelType core( dataPoints, options );
      
      // push the root into the queue
      working_list.push_back( std::make_pair( this,
                                              typename kernelType::State( &idx[0],
                                                                          static_cast<int>( idx.size() ),
                                                                          options ) ) );

      // initialize containers
      nodes.clear();

      // debugging:
      // int settledNum = 0;
      // int tmpCnt = 0;

      while ( !working_list.empty() ) {
        Tree<dataType,splitter> *node = fetch<order>(working_list).first;
        typename kernelType::State state = std::move( fetch<order>(working_list).second );
        pop<order>( working_list );
                  
        // emplace the current node
        node->nodeID = static_cast<int>( nodes.size() );
        nodes.push_back( NodeInfo() );
        nodes.back().node = node;
        
        // split
        std::vector<int> partition = std::move( core.split( state, node->judger ) );
        
        if ( 0 == partition[0] ) {
          // debugging:
          // DebugInfo( "%d -> %d + %d", state.len, partition[1], state.len - partition[1] );
          // ResumeOnRet();
          /* internal node */
          int branches = static_cast<int>( partition.size() )-1;
          node->child.resize(branches);
          for ( int i=0; i<branches; i++ ) {
            node->child[i].reset( new Tree() );
            working_list.push_back( std::make_pair( node->child[i].get(),
                                                    typename kernelType::State( state.idx + partition[i],
                                                                                partition[i+1] - partition[i],
                                                                                state ) ) );
          }
        } else {
          /* leaf node */
          for ( int i=0; i<state.len; i++ ) {
            nodes.back().store.push_back( state.idx[i] );
          }
          // debugging:
          // settledNum += state.len;
          // if ( settledNum > tmpCnt * 1000 ) {
          //   tmpCnt = settledNum / 1000 + 1;
          //   DebugInfo( "leaf: %d", settledNum );
          // }
        }
      }
    }

  private:
    /* ---------- private I/Os ---------- */
    void write( FILE* out ) const
    {
      judger.write( out );
      fwrite( &nodeID, sizeof(int), 1, out );
      
      unsigned char num = static_cast<unsigned char>( child.size() );
      fwrite( &num, sizeof(unsigned char), 1, out );
      for ( auto& every : child ) {
        every->write( out );
      }
    }

    Tree( FILE *in, std::vector<NodeInfo> &nodes )
    {
      judger.read( in );
      fread( &nodeID, sizeof(int), 1, in );
      nodes[nodeID].node = this;
      unsigned char num = 0;
      fread( &num, sizeof(unsigned char), 1, in );
      if ( 0 < num ) {
        child.resize( num );
        for ( int i=0; i<num; i++ ) {
          child[i].reset( new Tree( in, nodes ) );
        }
      }
    }

  public:
    /* ---------- public I/Os ---------- */
    void write( std::string filename ) const
    {
      WITH_OPEN( out, filename.c_str(), "w" );
      write( out );
      END_WITH( out );
    }

    Tree( std::string filename, std::vector<NodeInfo> &nodes )
    {
      WITH_OPEN( in, filename.c_str(), "r" );
      judger.read( in );
      fread( &nodeID, sizeof(int), 1, in );
      
      nodes[nodeID].node = this;
      unsigned char num = 0;
      fread( &num, sizeof(unsigned char), 1, in );
      if ( 0 < num ) {
        child.resize( num );
        for ( int i=0; i<num; i++ ) {
          child[i].reset( new Tree( in, nodes ) );
        }
      }
      END_WITH( in );
    }
    
    

    /* ---------- Properties ---------- */
    inline bool isLeaf() const
    {
      return 0 == child.size();
    }
    
    
    inline int levelSize( int level ) const
    {
      std::deque<Tree<dataType,splitter> > queue;
      if ( 0 == level || isLeaf() ) {
        return 1;
      }

      int s = 0;
      for ( auto& c : child ) {
        s += c->levelSize( level - 1 );
      }
      return s;
    }

    int depth() const
    {
      int max = 0;
      for ( auto& c : child ) {
        int t = c->depth();
        if ( t > max ) max = t;
      }
      return 1 + max;
    }

    const std::vector<std::unique_ptr< Tree<dataType,splitter> > >& getChild()
    {
      return child;
    }


    template <typename T>
    T reduce( std::function<T(const Tree<dataType,splitter>& node,const std::vector<T>&)> fun,
              std::function<T(const Tree<dataType,splitter>& node)> base ) const
    {
      if ( isLeaf() ) {
        return base( *this );
      } else {
        std::vector<T> results;
        for ( auto& c : child ) {
          results.push_back( c->reduce( fun, base ) );
        }
        return fun( *this, results );
      }
    }
                                    
              
              
              

    
    /* ---------- Queries ---------- */
    void collectLevel( int depth, std::vector<int>& store ) 
    {
      if ( 0 == depth ) {
        store.push_back( nodeID );
      } else {
        for ( auto& c : child ) {
          c->collectLevel( depth-1, store );
        }
      }
    }
    

    template <typename feature_t>
    int query( const feature_t &p ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );
      if ( isLeaf() ) {
        return nodeID;
      } else {
        return child[judger(p)]->query( p );
      }
    }
    
    template <typename feature_t>
    int query( const feature_t &p, int depth ) const
    {
      static_assert( std::is_same<typename ElementOf<feature_t>::type, dataType>::value,
                     "element of feature_t should have the same type as dataType." );            
      if ( depth == 0 || isLeaf() ) {
        return nodeID;
      } else {
        return child[judger(p)]->query( p, depth - 1 );
      }
    }

    /* ---------- Developer's Operations ---------- */
    const splitter<dataType> &getJudger() const
    {
      return judger;
    }
    
  };
}
