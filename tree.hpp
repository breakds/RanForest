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
    std::vector<kernel<dataType>::splitter> judger; // judgers of every node

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
    Forest() : dim(0), roots(), child(), judger(), store() {}

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
        judger.emplace_back();
        store.emplace_back();
      }

      while ( !worklist.empty() ) {
        size_t nodeId = fetch<order>(worklist).first;
        typename kernel<dataType>::State &state = fetch<order>(worklist).second;
        
        // try split
        judger[nodeId] = std::move( kernel<dataType>::ElectSplitter( dataPoints,
                                                                     dim,
                                                                     state ) );
        
                                    
      }
      

    }

                
               
              

               
    
               
    
    
  }
  
}
