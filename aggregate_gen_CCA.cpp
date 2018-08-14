//#include "Aggregate_Gen.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>
#include "aggregate_gen_SA.hpp"
using namespace Eigen;
using namespace std;


//return pos_sph matrix
MatrixXf aggregate_gen_CCA(const float& a, const int& num_sph_SA, const int& levels, const float& kf, const float& Df, const float& tol){
  int iter1max {1000};
  int iter2max {10000};

  int num_SA_agg= pow(2,levels);
  int totnum_sph= num_SA_agg*num_sph_SA;
  cout << "totnum_sph= " << totnum_sph << endl;

  vector<MatrixXf> pos_sph_SAagg;
  for(int i= 0; i<num_SA_agg; ++i ) pos_sph_SAagg.push_back(aggregate_gen_SA(a,num_sph_SA,kf,Df,tol));

  vector<MatrixXf> pos_sph_agg_old;

  // start time
  chrono::time_point<chrono::steady_clock> CCA_start_time= chrono::steady_clock::now();

  for(int k= 1; k <= levels; ++k){
      cout << "Generating aggregates in " << k << "th level ... please wait ..." << endl;
      vector<MatrixXf> pos_sph_agg;

      for(int L= 1; L <= pow(2,levels-k); ++L){
          //cout << "generating " << L << " th " << "aggregate in level" << k << endl;
          int N1= num_sph_SA*pow(2,k-1);
          int N2= N1;

          // aggregation process in k-th level
          MatrixXf pos_sph_agg1(3,N1);
          MatrixXf pos_sph_agg2(3,N2);

          if(k==1){
            pos_sph_agg1= pos_sph_SAagg[2*L-2];
            pos_sph_agg2= pos_sph_SAagg[2*L-1];
          } else {
            pos_sph_agg1= pos_sph_agg_old[2*L-2];
            pos_sph_agg2= pos_sph_agg_old[2*L-1];
          }

          // random rotation of aggregate 1 around the centroid (origin)
          float theta,phi; // polar angles
          theta= acos(2.0f*urf(re)-1.0f);
          phi= 2*M_PI*urf(re);

          Vector3f nvec; // unit vector of rotation axis
          nvec << sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta); //axis of rotation randomly choosen over unit sphere

          float psi; // angle of rotation
          psi= 2*M_PI*urf(re); // randomly choosen rotation angle

          Matrix3f rot_mat; // 3D rotation matrix
          rot_mat= AngleAxisf(psi,nvec);

          for(int i= 0; i < pos_sph_agg1.cols(); ++i) pos_sph_agg1.col(i)=rot_mat*pos_sph_agg1.col(i); // random rotation of agg1

          float R1_2; // square of gyration radius of agg1
          float R2_2; // square of gyration radius of agg2
          R1_2= a*a + pos_sph_agg1.colwise().squaredNorm().mean();
          R2_2= a*a + pos_sph_agg2.colwise().squaredNorm().mean();

          float rN;
          rN= a*a*(N1+N2)*(N1+N2)/(N1*N2)*pow((N1+N2)/kf,2/Df)-(N1+N2)/N2*R1_2-(N1+N2)/N1*R2_2;
          rN= sqrt(rN);

          bool found {false};

          //--- translation of agg1 to random position on the sphere of radius rN
          for(int iter1= 0; iter1 != iter1max; ++iter1){

              pos_sph_agg1.colwise() -= pos_sph_agg1.rowwise().mean(); // translate the centroid of agg1 to origin

              float phi;
              phi= 2*M_PI*urf(re);
              float u; // uniform random number [-1.0,1.0)
              u= 2*urf(re)-1.0f;

              Vector3f cen_pos_agg1; // centroid of agg1
              cen_pos_agg1 << sqrt(1.0f-u*u)*cos(phi), sqrt(1.0f-u*u)*sin(phi), u;
              cen_pos_agg1 *= rN; // cen_pos_agg1 is set to a uniform random cartesian position (x,y,z) over the spherical surface of radius rN
              pos_sph_agg1.colwise() += cen_pos_agg1;


              float theta;
              theta= acos(2.0f*urf(re)-1.0f);
              phi= 2*M_PI*urf(re);
              Vector3f nvec; // unit vector of rotation axis
              nvec << sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta); //axis of rotation randomly choosen over unit sphere




              // random rotation of aggregate 2 around the centroid (origin)
              #pragma omp parallel for
              for(int iter2= 0; iter2 < iter2max; ++iter2){
                  bool found_local;
                  #pragma omp atomic read
                  found_local= found;
                  if(found) continue;

                  MatrixXf pos_sph_agg1_thread(3,N1); // thread local copy of pos_sph_agg1
                  MatrixXf pos_sph_agg2_thread(3,N2); // thread local copy of pos_sph_agg2
                  pos_sph_agg1_thread= pos_sph_agg1;
                  pos_sph_agg2_thread= pos_sph_agg2;

                  float psi; // angle of rotation

                  #pragma omp critical
                  {
                     psi= 2*M_PI*urf(re);
                  }

                  Matrix3f rot_mat; // 3D rotation matrix
                  rot_mat= AngleAxisf(psi,nvec);
                  for(int i= 0; i < N2; ++i) pos_sph_agg2_thread.col(i)=rot_mat*pos_sph_agg2_thread.col(i); // random rotation of agg2
                  //  cout << "iter2 pos_sph_agg2 " << endl << pos_sph_agg2 << endl;

                  bool attached {false};
                  for(int i= 0; i < N2; ++i){
                      RowVectorXf dist_pair(N1);
                      dist_pair=(pos_sph_agg1_thread.colwise()-pos_sph_agg2_thread.col(i)).colwise().norm();
                    //  cout << "dist_pair" << endl << dist_pair << endl;
                      float mindist= dist_pair.minCoeff();
                      //cout << "2*a, mindist " << 2*a << ' ' << mindist << endl;
                      bool overrapped {mindist <= (2*a-tol)};
                      if(overrapped) {
                        //exit; // 20180813 removed
                        attached= false; // 20180813 added
                        break; // 20180813 added
                      } else {
                        attached= attached || (mindist < 2*a+tol);
                      }
                  }

                  // construction of new aggregate if attached condition is satisfied
                  if(attached){
                      //cout << "attached " << attached << endl;
                      MatrixXf pos_sph_agg_tmp(3,N1+N2);
                      pos_sph_agg_tmp.block(0,0,3,N1)= pos_sph_agg1_thread;
                      pos_sph_agg_tmp.block(0,N1,3,N2)= pos_sph_agg2_thread;
                      pos_sph_agg_tmp.colwise() -= pos_sph_agg_tmp.rowwise().mean(); // translate the centroid of new agg to origin
                      #pragma omp critical
                      {
                          pos_sph_agg.push_back(pos_sph_agg_tmp); // update the aggregate
                          found = true;
                          //break; // exit iter2 loop
                          if(k==levels)  cout << "number of threads for parallel search: " << omp_get_num_threads() << endl;
                      }

                  }
              }// iter2 loop

              if(found) break; // exit iter1 loop

              if(iter1 == iter1max){
                  cout << "calculation failed: iter1 reachs iter1max! " << endl;
                  abort();
              }
          }// iter1 loop
      }// L loop

      pos_sph_agg_old= move(pos_sph_agg);

  }// k loop

  auto CCA_end_time= chrono::steady_clock::now();

  cout << "CCA computation time: "
  << chrono::duration_cast<chrono::milliseconds>(CCA_end_time - CCA_start_time).count() << " milliseconds" << endl;

  MatrixXf pos_sph_agg_final;
  pos_sph_agg_final = pos_sph_agg_old[0];

  assert(totnum_sph == pos_sph_agg_final.cols());

  // final check of the attached condition
  double buf_factor= 5.0; //buffer factor of torelance
  bool non_overlapped;
  ArrayXf min_diff_from2a(totnum_sph);
  for(int i= 0; i< totnum_sph; ++i){
      ArrayXf distance_from_this_sph(totnum_sph);
      distance_from_this_sph= (pos_sph_agg_final.colwise()-pos_sph_agg_final.col(i)).colwise().norm().array();
      non_overlapped= {(distance_from_this_sph-(2*a-buf_factor*tol) > 0)(i) == false};
      for(int j= 0; j< totnum_sph; ++j){
          if(j != i) {
            non_overlapped= non_overlapped && ((distance_from_this_sph-(2*a-buf_factor*tol) > 0)(j) == true);
          }
      }
     //cout << i << ", " << (distance_from_this_sph-(2*a-buf_factor*tol) >= 0).transpose() << endl ;
      min_diff_from2a(i)=((pos_sph_agg_final.colwise()-pos_sph_agg_final.col(i)).colwise().norm().array()-2*a).array().abs().minCoeff();
  }
  bool attached_at_least_1pair {min_diff_from2a.maxCoeff() < buf_factor*tol};

  string non_overlapped_test= {non_overlapped ? "Success" : "Failed"};
  string attached_at_least_1pair_test= {attached_at_least_1pair ? "Success" : "Failed"};
  cout << "non_overlapped_check: " << non_overlapped_test << endl;
  cout << "attached_at_least_1pair_check: " << attached_at_least_1pair_test << endl;
  if(non_overlapped && attached_at_least_1pair){
    return pos_sph_agg_final;
  }else{
    cout << "Geometry check failed! Please try other (kf, Df) pair." << endl;
    exit(1);
  }



}
