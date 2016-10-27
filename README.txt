**************************************************************************************
This software is released under the MIT License, see LICENSE.txt
**************************************************************************************

Code title:
  "aggregate_gen" -- An efficient C++ code for generating fractal cluster of spheres


General Descriptions:
  A C++ code for calculating the geometrical coordinate of fractal-like cluster of spheres, using the tunable cluster-cluster aggregation (CCA) method.
  A user is allowed to change the geometrical parameters: the number of spherules N, the fractal prefactor kf, and the fractal dimension Df.


Developer's info:
  Nobuhiro Moteki,
  Assistant Professor,Department of Earth and Planetary Sciences, The University of Tokyo
  moteki@eps.s.u-tokyo.ac.jp


Update history:
  version 1.0.1 (October 27, 2016);  Initial release


Prerequisites:
  1. A C++ compiler supporting C++11. The author tested the code using the gcc version 6.2.0.
  2. A C++ library for linear algebra "Eigen" with version 3.2 or newer. The author tested the code using the Eigen version 3.2.10.


Usage:
  1. Modify the compilation parameters in "Makefile" according to your environment. In particular,
      please change the "CXX" (= C++ complier) and "INCLUDES" (= absolute path of Eigen folder).

  2. Compilation and linking are performed by a command
      > make

  3. Run the executable by a command
      > ./aggregate_gen [parameter_list]

      Here, the [parameter_list] denotes a set of space-delimited input parameters as follows.
      [parameter_list] == num_sph_SA levels kf Df
          num_sph_SA: The number of spherules in initial clusters in the 1st level of cluster-cluster aggregation (4).
          levels: The number of hierarchical levels in cluster-cluster aggregation (10).
          kf: The fractal prefactor (1.3).
          Df: The fractal dimension (2.5).
      If you skip the [parameter_list], the code uses the default values in the parenthesis.
      Total number of spherules is num_sph_SA*(2^levels).

  4. Output data
      An execution of aggregate_gen outputs a text file named "agg_N***_kf***_Df***.out", where the three strings
      *** indicate the total number of spherules N, the fractal prefactor kf, and the fractal dimension Df.
      Each row shows the cartesian coordinate "x y z" of each spherules relative to the centroid of cluster, where the
      length is scaled such that spherule radius is 1.0.


Other Notes:
  By default, the seed of random_number sequence is automatically changed in each execution (using random_device{}()).
  If necessary, the user can fix the seed by commenting the line 14 and uncommenting line 15 in "aggregate_gen_main.cpp".


Theoretical basis:
  This code uses the mathematical theory described in "A. V. Filippov, M. Zurita, and D. E. Rosner , Journal of Colloid and Interface Science 229, 261–273 (2000)".
