#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "aggregate_gen_SA.hpp"
using namespace Eigen;
using namespace std;

//return pos_sph matrix
MatrixXf aggregate_gen_SA(const float& a, const int& num_sph, const float& kf, const float& Df, const float& tol){
    MatrixXf pos_sph= MatrixXf::Zero(3,num_sph);
    pos_sph.col(0) << -a,0,0;
    pos_sph.col(1) << a,0,0;

    Vector3f pos_cand(3);
    float rN;
    float phi;
    float u;
    ArrayXf dist_sph;

    for(int N=3; N <= num_sph; ++N){
        rN= N*N*a*a/(N-1)*pow(N/kf,2/Df)-N*a*a/(N-1)-N*a*a*pow((N-1)/kf,2/Df);
        rN= sqrt(rN);
        dist_sph.resize(N-1);

        while(1){
            phi= 2*M_PI*urf(re);
            u= 2*urf(re)-1;
            pos_cand << sqrt(1.0-u*u)*cos(phi), sqrt(1.0-u*u)*sin(phi), u;
            pos_cand *= rN;

            for(int i= 0; i < N-1; ++i){
                dist_sph(i)=(pos_sph.col(i)-pos_cand).norm();
            }

            if(dist_sph.minCoeff() < 2*a+tol && dist_sph.minCoeff() > 2*a-tol){
                pos_sph.col(N-1)= pos_cand;
                pos_sph.colwise() -= pos_sph.rowwise().mean(); // translate the centroid to origin
                break;
            }
        }
    }

    // final check of the attached condition
    ArrayXf mindist_c(num_sph);
    for(int i= 0; i< num_sph; ++i){
        mindist_c(i)=((pos_sph.colwise()-pos_sph.col(i)).colwise().norm().array()-2*a).array().abs().minCoeff();
    }
    assert(mindist_c.maxCoeff() < tol);

    return pos_sph;
}
