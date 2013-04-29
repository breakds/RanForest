#pragma once

#include <type_traits>


namespace ran_forest
{
  template <typename dataType = float, template <typename> class kernel>
  class Forest
  {
  private:
    int dim; // dimension of feature vectors
    std::vector<size_t> roots; // nodeID of roots
    std::vector<std::vector<size_t> > child; // children IDs for every node
    std::vector<kernel<dataType>::splitter> judge; // judges of every node

    // training datapoints IDs for every node. currently only valid
    // for *leaf* nodes
    std::vector<std::vector<size_t> > store; 



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
    
    // constructor 0: default constructor
    Forest() : dim(0), roots(), child(), judge(), store() {}

    // constructor 1: reading constructor


    // grow the forest


    // grow one tree, return the nodeID of the root
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
#pragma omp critical
      {
        root = child.size();
        child.emplace_back();
        judge.emplace_back();
        store.emplace_back();
      }

      worklist.push_back( std::make_pair( root,
                                          typename kernel<dataType>::State( &idx[0],
                                                                            idx.size(),
                                                                            options ) ) );
      std::vector<int> label( dataPoints.size() );


      while ( !worklist.empty() ) {
        size_t nodeId = fetch<order>(worklist).first;
        typename kernel<dataType>::State state = fetch<order>(worklist).second;
        pop<order>( worklist );
        
        // try split
        ElectionStatus status = kernel<dataType>::ElectSplitter( dataPoints, dim, state, judge[nodeID] );
        int maxLabel = -1;
        if ( SUCCESS == status ) {
          // calculate branch label
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
            int i = curpos[k];
            while ( i < partition[k+1] ) {
              int k1 = label[state.idx[i]];
              if ( k1 != k ) {
                j = curpos[k1]++;
                size_t tmp = state.idx[i];
                state.idx[i] = state.idx[j];
                state.idx[j] = tmp;
              } else {
                i++;
              }
            }
          }

          // split
#pragma   omp critical
          {
            for ( int k=0; k<=maxLabel; k++ ) {
              size_t id = child.size();
              child.emplace_back();
              judge.emplace_back();
              store.emplace_back();
              child[nodeID].push_back( id );
              worklist.push_back( std::make_pair( id,
                                                  typename kernelType::State( state.idx + partition[k],
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
      

    }

                
               
              

               
    
               
    
    
  }
  
}
