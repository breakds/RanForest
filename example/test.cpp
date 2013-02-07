#include <cstdio>
#include <random>
#include <vector>
#include "../component/forest.hpp"

using namespace ran_forest;

class multi_normal_gen
{
  std::default_random_engine engine;
  std::vector<std::normal_distribution<float> > generators;
public:
  template <typename meanArrayType>
  multi_normal_gen( const meanArrayType& means, float stddev, int dim )
  {
    for ( int i=0; i<dim; i++ ) {
      generators.emplace( generators.end(), means[i], stddev );
    }
  }
  
  void operator()( float* v )
  {
    int i = 0;
    for ( auto& gen : generators ) {
      v[i++] = gen(engine);
    }
  }
};



int main()
{
  int K = 3;
  int perClass = 500;
  int dim = 50;
  float stddev = 0.05;
  std::vector<std::vector<float> > centers(K);


  // Generate centers
  std::default_random_engine engine;
  std::uniform_real_distribution<float> dist( 0.0, 1.0 );
  for ( int k=0; k<K; k++ ) {
    for ( int i=0; i<dim; i++ ) {
      centers[k].push_back( dist(engine) );
    }
  }

  std::vector<std::vector<float> > features;
  for ( int k=0; k<K; k++ ) {
    multi_normal_gen gen( centers[k], stddev, dim );
    for ( int i=0; i<perClass; i++ ) {
      features.emplace( features.end() );
      features.back().resize(dim);
      gen( &features.back()[0] );
    }
  }

  std::vector<int> idx(K*perClass);
  for ( int i=0; i<K*perClass; i++ ) {
    idx[i] = i;
  }

  typename SimpleKernel<std::vector<float>, BinaryOnAxis>::Options options;

  
  options.dim = dim;
  options.converge = 0.1;
  options.stopNum = 10;


  // Tree<float,BinaryOnAxis> tree;
  // tree.grow<SimpleKernel>( features, idx, nodes, options );

  Forest<float,BinaryOnAxis> forest;
  forest.grow<SimpleKernel>( 10, features, options, 0.5 );


  int leaves = 0;
  int count = 0;
  for ( auto& node : forest ) {
    if ( 0 < node.store.size() ) {
      leaves++;
      int c = node.store[0] / perClass;
      for ( auto& ele : node.store ) {
        if ( c != (ele / perClass) ) {
          count--;
          break;
        }
      }
      count++;
    }
  }
  Info( "%d/%d pass", count, leaves );
  

  count = 0;
  for ( int i=0; i<K*perClass; i++ ) {
    std::vector<int> nodeIDs = forest.query( features[i] );
    for ( int& nodeID : nodeIDs ) {
      for ( auto& ele : forest[nodeID].store ) {
        if ( ele == i ) {
          count++;
          break;
        }
      }
    }
  }
  Info( "%d/%d pass", count, K * perClass );

  
  
  return 0;
}
