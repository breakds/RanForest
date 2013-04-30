#include <cstdio>
#include <random>
#include <vector>
#include "../RanForest.hpp"

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

  
  VP<float>::Options options;
  options.converge = 0.01;
  

  Forest<float,VP> forest;
  forest.grow( 1, features, dim, options );


  int count = 0;
  for ( int i=0; i<K*perClass; i++ ) {
    std::vector<size_t> nodeIDs = forest.query( features[i] );
    for ( size_t& nodeID : nodeIDs ) {
      for ( auto& ele : forest.getStore( nodeID ) ) {
        if ( ele == static_cast<size_t>( i ) ) {
          count++;
          break;
        }
      }
    }
  }
  Info( "%d/%d pass", count, K * perClass );

  
  
  return 0;
}
