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
vector<MatrixXd> aggregate_gen_CCA(const double& a, const int& num_sph_SA, const int& levels, const double& kf, const double& Df, const double& tol){

  vector<MatrixXd> pos_sph_agg_output; // vector of aggregate geometry from level 0 to "levels".
  int iter1max {1000};
  int iter2max {10000};

  int num_SA_agg= pow(2,levels);
  int final_totnum_sph= num_SA_agg*num_sph_SA;
  cout << "final totnum_sph= " << final_totnum_sph << endl;

  vector<MatrixXd> pos_sph_SAagg;
  for(int i= 0; i<num_SA_agg; ++i ) pos_sph_SAagg.push_back(aggregate_gen_SA(a,num_sph_SA,kf,Df,tol));

  pos_sph_agg_output.push_back(pos_sph_SAagg[0]); // level 0 aggregate


  vector<MatrixXd> pos_sph_agg_old;

  // start time
  chrono::time_point<chrono::steady_clock> CCA_start_time= chrono::steady_clock::now();

  for(int k= 1; k <= levels; ++k){
      cout << "Generating aggregates in " << k << "th level ... please wait ..." << endl;
      vector<MatrixXd> pos_sph_agg;

      for(int L= 1; L <= pow(2,levels-k); ++L){
          //cout << "generating " << L << " th " << "aggregate in level" << k << endl;
          int N1= num_sph_SA*pow(2,k-1);
          int N2= N1;

          // aggregation process in k-th level
          MatrixXd pos_sph_agg1(3,N1);
          MatrixXd pos_sph_agg2(3,N2);

          if(k==1){
            pos_sph_agg1= pos_sph_SAagg[2*L-2];
            pos_sph_agg2= pos_sph_SAagg[2*L-1];
          } else {
            pos_sph_agg1= pos_sph_agg_old[2*L-2];
            pos_sph_agg2= pos_sph_agg_old[2*L-1];
          }

          // random rotation of aggregate 1 around the centroid (origin)
          double theta,phi; // polar angles
          theta= acos(2.0f*urf(re)-1.0f);
          phi= 2*M_PI*urf(re);

          Vector3d nvec; // unit vector of rotation axis
          nvec << sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta); //axis of rotation randomly choosen over unit sphere

          double psi; // angle of rotation
          psi= 2*M_PI*urf(re); // randomly choosen rotation angle

          Matrix3d rot_mat; // 3D rotation matrix
          rot_mat= AngleAxisd(psi,nvec);

          for(int i= 0; i < pos_sph_agg1.cols(); ++i) pos_sph_agg1.col(i)=rot_mat*pos_sph_agg1.col(i); // random rotation of agg1

          double R1_2; // square of gyration radius of agg1
          double R2_2; // square of gyration radius of agg2
          R1_2= a*a + pos_sph_agg1.colwise().squaredNorm().mean();
          R2_2= a*a + pos_sph_agg2.colwise().squaredNorm().mean();

          double rN;
          rN= a*a*(N1+N2)*(N1+N2)/(N1*N2)*pow((N1+N2)/kf,2/Df)-(N1+N2)/N2*R1_2-(N1+N2)/N1*R2_2;
          rN= sqrt(rN);

          bool found {false};

          //--- translation of agg1 to random position on the sphere of radius rN
          for(int iter1= 0; iter1 != iter1max; ++iter1){

              pos_sph_agg1.colwise() -= pos_sph_agg1.rowwise().mean(); // translate the centroid of agg1 to origin

              double phi;
              phi= 2*M_PI*urf(re);
              double u; // uniform random number [-1.0,1.0)
              u= 2*urf(re)-1.0;

              Vector3d cen_pos_agg1; // centroid of agg1
              cen_pos_agg1 << sqrt(1.0-u*u)*cos(phi), sqrt(1.0-u*u)*sin(phi), u;
              cen_pos_agg1 *= rN; // cen_pos_agg1 is set to a uniform random cartesian position (x,y,z) over the spherical surface of radius rN
              pos_sph_agg1.colwise() += cen_pos_agg1;

              double theta;
              theta= acos(2.0*urf(re)-1.0);
              phi= 2*M_PI*urf(re);
              Vector3d nvec; // unit vector of rotation axis
              nvec << sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta); //axis of rotation randomly choosen over unit sphere

              // random rotation of aggregate 2 around the centroid (origin)
              #pragma omp parallel for
              for(int iter2= 0; iter2 < iter2max; ++iter2){
                  bool found_local;
                  #pragma omp atomic read
                  found_local= found;
                  if(found) continue;

                  MatrixXd pos_sph_agg1_thread(3,N1); // thread local copy of pos_sph_agg1
                  MatrixXd pos_sph_agg2_thread(3,N2); // thread local copy of pos_sph_agg2
                  pos_sph_agg1_thread= pos_sph_agg1;
                  pos_sph_agg2_thread= pos_sph_agg2;

                  double psi; // angle of rotation

                  #pragma omp critical
                  {
                     psi= 2*M_PI*urf(re);
                  }

                  Matrix3d rot_mat; // 3D rotation matrix
                  rot_mat= AngleAxisd(psi,nvec);
                  for(int i= 0; i < N2; ++i) pos_sph_agg2_thread.col(i)=rot_mat*pos_sph_agg2_thread.col(i); // random rotation of agg2
                  //  cout << "iter2 pos_sph_agg2 " << endl << pos_sph_agg2 << endl;

                  bool attached {false};
                  for(int i= 0; i < N2; ++i){
                      RowVectorXd dist_pair(N1);
                      dist_pair=(pos_sph_agg1_thread.colwise()-pos_sph_agg2_thread.col(i)).colwise().norm();
                    //  cout << "dist_pair" << endl << dist_pair << endl;
                      double mindist= dist_pair.minCoeff();
                      //cout << "2*a, mindist " << 2*a << ' ' << mindist << endl;
                      bool overrapped {mindist <= 2*a};
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
                      MatrixXd pos_sph_agg_tmp(3,N1+N2);
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
      pos_sph_agg_output.push_back(pos_sph_agg_old[0]);

  }// k loop

  auto CCA_end_time= chrono::steady_clock::now();

  cout << "CCA computation time: "
  << chrono::duration_cast<chrono::milliseconds>(CCA_end_time - CCA_start_time).count() << " milliseconds" << endl;

  // geometry check for output
  bool non_overlapped_all {true};
  bool attached_at_least_1pair_all {true};

  for(int k= 0; k <= levels; ++k){
      MatrixXd pos_sph_agg_k_out {pos_sph_agg_output[k]};
      int num_sph = num_sph_SA*pow(2,k);
      assert(num_sph == pos_sph_agg_k_out.cols());
      // final check of the attached condition
      bool non_overlapped;
      ArrayXd min_diff_from2a(num_sph);
      for(int i= 0; i< num_sph; ++i){
          ArrayXd distance_from_this_sph(num_sph);
          distance_from_this_sph= (pos_sph_agg_k_out.colwise()-pos_sph_agg_k_out.col(i)).colwise().norm().array();
          non_overlapped= {(distance_from_this_sph-2*a >= 0)(i) == false};
          for(int j= 0; j< num_sph; ++j){
              if(j != i) {
                non_overlapped= non_overlapped && ((distance_from_this_sph-2*a > 0)(j) == true);
              }
          }
         //cout << i << ", " << (distance_from_this_sph-(2*a-buf_factor*tol) >= 0).transpose() << endl ;
          min_diff_from2a(i)=((pos_sph_agg_k_out.colwise()-pos_sph_agg_k_out.col(i)).colwise().norm().array()-2*a).array().abs().minCoeff();
      }
      bool attached_at_least_1pair {min_diff_from2a.maxCoeff() < tol};
      string non_overlapped_test= {non_overlapped ? "Success" : "Failed"};
      string attached_at_least_1pair_test= {attached_at_least_1pair ? "Success" : "Failed"};
      cout << k << "th level, non_overlapped_check         : " << non_overlapped_test << endl;
      cout << k << "th level, attached_at_least_1pair_check: " << attached_at_least_1pair_test << endl;

      non_overlapped_all = non_overlapped_all && non_overlapped;
      attached_at_least_1pair_all = attached_at_least_1pair_all && attached_at_least_1pair;
  }

  if(non_overlapped_all && attached_at_least_1pair_all){
    return pos_sph_agg_output;
  }else{
    cout << "Geometry check failed! Please change tol." << endl;
    exit(1);
  }



}
