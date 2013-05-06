// This file is part of RanForest, a lightweight C++ template library
// for random forest.
//
// By BreakDS <breakds@cs.wisc.edu> - http://www.unlicense.org/ (public domain)
// 
//
//
// This file implements the quasi K-Means clustering algorithm, as
// described below:
//
// 1. The input of the algorithm is a preliminary clustering result,
//    namely a bipartite graph, where each data point is associated
//    with one or several clusters.
//    
// 2. The output of the algorithm is also a bipartite graph.
//
// 3. On every iteration, each data point will be associated to the
//    top K nearest centers (clusters) among all the associated
//    centers in the original graph, and each center (cluster) will be
//    updated as the mean of all associated data points.




#pragma once

#include <vector>
#include <memory>
#include "../RanForest.hpp"
#include "LLPack/algorithms/heap.hpp"


using ran_forest::Bipartite;

namespace ran_forest {

  template <typename dataType = float>
  class TMeanShell
  {
  public:

    // the dimension of the data points that will be handled
    int dim;

    // the centers of the resulting clusters
    std::vector<std::vector<dataType> > centers;

    struct Options
    {
      // the maximum number of iterations
      int maxIter;

      // the maximum number of clusters that one data point can
      // connect to
      int replicate;

      // if the energy function change is below this level, the
      // algorithm stops
      double converge;

      // the bandwidth used to calculate the exponential
      // transformation on distances between data points and clusters,
      // whose normalized version will be used as the edge weight of
      // the resulting bipartite graph
      double wtBandwidth;
      
      Options() : maxIter(20), replicate(10), converge(1e-5), wtBandwidth(100.0) {}
    } options;

    TMeanShell( int dimension ) : dim(dimension), centers(), options() {}


    // +-------------------------------------------------------------------------------
    // Input and Output
  private:

    // Put a sign at the end of a binary file
    static void seal( FILE* out )
    {
      char ch[10] = "ENDTMEANS";
      fwrite( ch, sizeof(char), 10, out );
    }
    
    // Check whether a correct sign is present 
    static bool unseal( FILE* in )
    {
      char ch[10];
      fread( ch, sizeof(char), 10, in );
      return 0 == strcmp( ch, "ENDTMEANS" );
    }


  public:

    
    void write( std::string filename )
    {
      WITH_OPEN( out, filename.c_str(), "wb" );
      size_t len = centers.size();
      fwrite( &len, sizeof(size_t), 1, out );
      for ( auto& ele : centers ) {
        fwrite( &ele[0], sizeof(dataType), dim, out );
      }
      seal( out );
      END_WITH( out );
    }


    void read( std::string filename )
    {
      WITH_OPEN( in, filename.c_str(), "rb" );
      size_t len = 0;
      fread( &len, sizeof(size_t), 1, in );
      centers.resize( len );
      for ( auto& ele : centers ) {
        ele.size( dim );
        fread( &ele[0], sizeof(dataType), dim, in );
      }
      if ( !unseal( in ) ) {
        Error( "TMeanShell: unseal() failed, might be due to wrong centers data." );
        exit( -1 );
      }
      END_WITH( in );
    }

  private:

    template <typename feature_t, template <typename T = feature_t, typename... restArgs> class container>
    void CenterMeans( std::vector<std::vector<dataType> > &centers,
                      const container<feature_t>& feat,
                      Bipartite& n_to_l )
    {
      
      size_t L = n_to_l.sizeB();

#     pragma omp parallel for
      for ( size_t l=0; l<L; l++ ) {
        std::fill( centers[l].begin(), centers[l].end(), 0.0 );
        auto& _to_n = n_to_l.to( l );
        if ( 0 < _to_n.size() ) {
          size_t count = 0;
          for ( auto& ele : _to_n ) {
            count++;
            size_t n = ele.first;
            for ( int j=0; j<dim; j++ ) {
              centers[l][j] += feat[n][j];
            }
          }
          double wt = static_cast<float>( 1.0 / count );
          for ( auto& component : centers[l] ) {
            component *= wt;
          }
        }
      }
    }
    


  public:

  

    template <typename feature_t, template <typename T = feature_t, typename... restArgs> class container>
    void Clustering( const container<feature_t> &feat,
                     Bipartite& n_to_l,
                     bool silent = false )
    {
      size_t N = n_to_l.sizeA();
      size_t L = n_to_l.sizeB();
      
      centers.resize( L );
      for ( auto& ele : centers ) {
        ele.resize( dim );
      }

      CenterMeans( centers, feat, n_to_l );

      Bipartite bimap( N, L );

      double lastEnergy = 0.0;
              
      for ( int iter=0; iter<options.maxIter; iter++ ) {
        

        bimap.clear();

        if ( !silent ) {
          Info( "TMeans iter %d", iter );
        }
        // pick centers
#       pragma omp parallel for
        for ( size_t n=0; n<N; n++ ) {
          auto& _to_l = n_to_l.from( n );
          heap<double,size_t> ranker( options.replicate );
          for ( auto& ele : _to_l ) {
            size_t l = ele.first;
            double dist = 0.0;
            for ( int j=0; j<dim; j++ ) {
              double tmp = centers[l][j] - feat[n][j];
              dist += tmp * tmp;
            }
            ranker.add( dist, l );
          }
          
#         pragma omp critical
          for ( int j=0; j<ranker.len; j++ ) {
            bimap.add( n, ranker[j], 1.0 / options.replicate );
          }
          
        } // end for n

        CenterMeans( centers, feat, bimap );



        // Calculate Energy
        double energy = 0.0;
#       pragma omp parallel for reduction(+ : energy)
        for ( size_t n=0; n<N; n++ ) {
          auto& _to_l = bimap.from( n );
          for ( auto& ele : _to_l ) {
            size_t l = ele.first;
            for ( int j=0; j<dim; j++ ) {
              double tmp = centers[l][j] - feat[n][j];
              energy += tmp * tmp;
            }
          }
        }
        
        
        if ( 0 <iter && fabs(lastEnergy-energy) < options.converge ) {
          break;
        }
        
        lastEnergy = energy;

        if ( !silent ) {
          printf( "Energy: %.5lf\n", energy );
        }
        
      } // end for iter

      // update alphas
#     pragma omp parallel for
      for ( size_t n=0; n<N; n++ ) {
        auto& _to_l = bimap.getSetFrom( n );
        if ( 0 < _to_l.size() ) {

          double s = 0.0;
          for ( auto& ele : _to_l ) {
            size_t l = ele.first;

            // calculate l2 distance
            double dist = 0.0;
            for ( int j=0; j<dim; j++ ) {
              double tmp = centers[l][j] - feat[n][j];
              dist += tmp * tmp;
            }
            dist = sqrt( dist );

            ele.second = exp( - dist / options.wtBandwidth );
	    
            s += ele.second;

          }
          
          s = 1.0 / s;
          for ( auto& ele : _to_l ) {
            ele.second *= s;
          }
        }
      }
      n_to_l = std::move( bimap );
    }


    template <typename feature_t>
    inline void concentrate( const feature_t &p, std::vector<std::pair<int,double> > &membership ) const
    {
      heap<double,size_t> ranker( options.replicate );
      for ( auto& ele : membership ) {
        size_t l = ele.first;
        double dist = 0.0;
        for ( int j=0; j<options.dim; j++ ) {
          double tmp = centers[l][j] - p[j];
          dist += tmp * tmp;
        }
        ranker.add( sqrt( dist ), l );
      }
      membership.resize( ranker.len );
      double s = 0.0;
      for ( int i=0; i<ranker.len; i++ ) {
        membership[i].first = ranker[i];
        membership[i].second = exp( - ranker(i) / options.wtBandwidth );
	s += membership[i].second;
      }
      // normalization
      s = 1.0 / s;
      for ( int i=0; i<ranker.len; i++ ) {
	membership[i].second *= s;
      }
    }
  };
}
