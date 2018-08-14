
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <string>
#include <Eigen/Dense>
#include "aggregate_gen_CCA.hpp"
using namespace Eigen;
using namespace std;

default_random_engine re(random_device{}()); // random seed
//default_random_engine re; //default fixed seed
uniform_real_distribution<float> urf {0,1.0f};


int main(int argc, char *argv[]){

    // argv[0] == program name
    // argv[1] == num_sph_SA
    // argv[2] == levels
    // argv[3] == fractal prefactor kf
    // argv[4] == fractal dimension Df
    // argv[5] == error tolerance of attached distance

    float a= {1.0}; // radius of spherule
    float tol= {1.0e-3}; // error torelance of the distance between attached spherules
    int num_sph_SA;
    int levels;
    float kf;
    float Df;

    switch (argc)
    {
    case 1:
        // use default parameter_list
        num_sph_SA= {4}; // number of spherule in aggregate generated given by SA
        levels= {10}; //number of hierarchical coagulation
        kf= {1.0}; // fractal prefactor
        Df= {2.5}; // fractal dimension
        tol= {1.0e-3}; // error torelance of the distance between attached spherules
        break;
    case 5:
        // parameter_list from the command line arguments
        num_sph_SA= atoi(argv[1]); // number of spherule in aggregate generated given by SA
        levels= atoi(argv[2]); //number of hierarchical coagulation
        kf= atof(argv[3]); // fractal prefactor
        Df= atof(argv[4]); // fractal dimension
        break;
    default:
        runtime_error("incorrect parameter_list");
        break;
    }

    cout << "num_sph_SA: " << num_sph_SA << endl;
    cout << "levels: " << levels << endl;
    cout << "kf: " << kf << endl;
    cout << "Df: " << Df << endl;
    cout << "tol: " << tol << endl;

    MatrixXf pos_sph= aggregate_gen_CCA(a,num_sph_SA,levels,kf,Df,tol);
    int totnum_sph= pos_sph.cols();

    cout << "totnum_sph: " << totnum_sph << endl;

    string oname;
    oname = "agg_N"+to_string(totnum_sph)+"_kf"+to_string(kf).erase(3)+"_Df"+to_string(Df).erase(3)+".out";
    ofstream ost {oname};
    for(int i=0; i<totnum_sph; ++i){
        ost << scientific << showpoint
         //<< a << ' '
         << pos_sph(0,i) << ' '
         << pos_sph(1,i) << ' '
         << pos_sph(2,i) << endl;
    }

}
