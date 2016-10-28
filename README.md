##  "aggregate_gen" -- An efficient C++ code for generating fractal cluster of spheres
---

### General Descriptions:
  A C++ code for calculating the geometrical coordinate of fractal-like cluster of spheres, using the tunable cluster-cluster aggregation (CCA) method.
  A user is allowed to change the geometrical parameters: the number of spherules N, the fractal prefactor kf, and the fractal dimension Df.


### Prerequisites:
  1. A C++ compiler supporting C++11. The author tested the code using the gcc version 6.2.0.
  2. A C++ library for linear algebra "Eigen" with version 3.2 or newer. The author tested the code using the Eigen version 3.2.10.


### Usage:
  1. Modify the compilation parameters in "Makefile" according to your environment.
  In particular,please change the "CXX" (= C++ complier) and "INCLUDES" (= absolute path of Eigen folder).

  2. Compilation and linking are performed by

            make

  3. Run the executable by

           ./aggregate_gen [parameter_list]

    [parameter_list] == *num_sph_SA* *levels*  *kf*  *Df*

    * *num_sph_SA* : Number of spherules in initial clusters in the 1st level of cluster-cluster aggregation (4).
    * *levels* : Number of hierarchical levels in cluster-cluster aggregation (10).
    * *kf* : Fractal prefactor (1.3).
    * *Df* : Fractal dimension (2.5).

    If you skip the [parameter_list], the code uses the default values in the parenthesis. Total number of spherules is *N* = num_sph_SA*2<sup>levels</sup>.


### Output data
  The result is written in a text file "agg_N???_kf???_Df???.out", where the three strings ??? indicate the total number of spherules *N*, the fractal prefactor *kf*, and the fractal dimension *Df*.

  In this file, each row shows the cartesian coordinate "x y z" of each spherules relative to the centroid of cluster, where the
  length is scaled such that spherule radius is 1.0.

### Accuracy
Maximum error of the center-to-center distance between two-spherules in contact never exceeds 0.005 (typically < 0.001).

### Performance
This code is probably one of the fastest tunable-CCA codes reported so far. For computing a fractal-like aggregate with parameter Df=2.1, kf=1.1, and N=65536, this code consumed only 3784 seconds (Intel Xeon E5 (3.5GHz)). In comparison, Skorupski et al. 2014<sup>(a)</sup> reported that their original tunable-CCA code consumed 21809 seconds (AMD athron II X4 640 (3.0GHz)) for computing an aggregate with the same parameter.

<sup>(a)</sup> Skorupski et al. 2014, A fast and accurate implementation of tunable algorithms used for generation of fractal-like aggregate models, *Physica A* 404, 106-117.

### Random number seed
If necessary, random number seed can be changed at line14-15 of "aggregate_gen_main.cpp".

    (line 14) default_random_engine re(random_device{}()); // random seed
    (line 15) //default_random_engine re; // fixed seed

By default, the seed of random_number sequence is automatically changed in each execution (using random_device{}()). Fixed seed will be used if you uncomment 15th line and comment out the 14th line.

### Theory
  This code uses the mathematical theory described in "Filippov, A.V. et al. 2000, *Journal of Colloid and Interface Science* 229, 261–273".

### Contact
nobuhiro.moteki@gmail.com

### Licence
Please see LICENSE.txt
