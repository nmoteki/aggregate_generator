##  "aggregate_gen" -- An efficient C++ code for generating fractal cluster of spheres
---

### History:
#### Aug 14, 2018:
  * Bugfix: An algorithm error causing the partial overlapping between some spherules has been fixed.
  * Improve reliability: Internal tests of non-overlapping and attached conditions are added.
  * Improve performance: The MonteCarlo search of an attached-and-nonoverlapped aggregate configuration is parallelized using openmp.

#### Oct 28, 2016:
  * Initial release.

### General Descriptions:
  * A C++ code for calculating the geometrical coordinate of fractal-like cluster of spheres, using the tunable cluster-cluster aggregation (CCA) method.
  * A user is allowed to change the geometrical parameters: the number of spherules N, the fractal prefactor kf, and the fractal dimension Df.
  * The excecusion speed automatically improves with the available number of CPU cores.


### Prerequisites:
  1. A C++ compiler supporting C++17 and openmp. The author tested the code using the g++ version 7.3.0.
  2. A C++ library for linear algebra "Eigen" with version 3.3.5 or newer. The author tested the code using the Eigen version 3.3.5.
  3. An automated compilation and build system "CMake". The author tested the code using the CMake version 3.10.2.


### Usage:
  1. Modify the parameters in "CMakeLists.txt" according to your environment.
  In particular,please change the "CMAKE_CXX_COMPILER" (= C++ complier) and "include_directories" (= absolute path of Eigen folder).

  2. Compilation and linking are performed by
```
            cmake .
            make
```

  3. Run the executable by
```
           ./aggregate_gen_main [parameter_list]
```
    [parameter_list] == *num_sph_SA* *levels*  *kf*  *Df*

    * *num_sph_SA* : Number of spherules in initial clusters in the 1st level of cluster-cluster aggregation (4).
    * *levels* : Number of hierarchical levels in cluster-cluster aggregation (10).
    * *kf* : Fractal prefactor (1.0).
    * *Df* : Fractal dimension (2.5).

    If you skip the [parameter_list], the code uses the default values in the parenthesis. Total number of spherules is *N* = num_sph_SA*2<sup>levels</sup>.


### Output data
  The result is written in a text file "agg_N???_kf???_Df???.out", where the three strings ??? indicate the total number of spherules *N*, the fractal prefactor *kf*, and the fractal dimension *Df*.

  In this file, each row shows the cartesian coordinate "x y z" of each spherules relative to the centroid of cluster, where the
  length is scaled such that spherule radius is 1.0.

### Accuracy
Maximum error of the center-to-center distance between two-spherules in contact never exceeds 0.005 (typically < 0.001). You can change the tolerance by changing the parameter "tol" defined in "aggregate_gen_main.cpp"

### Limitation of (kf,Df)
Under the constraint of the fractal-scaling law, finding of an attached-and-nonoverlapped aggregate configuration could be unrealistic for compact (i.e., closely-packed) aggregates. To avoid the infinite looping of the search algorithm, please choose the geometrical parameters (kf,Df) to be kf+Df < ~3.4.

### Random number seed
If necessary, random number seed can be changed at line14-15 of "aggregate_gen_main.cpp".

    (line 14) default_random_engine re(random_device{}()); // random seed
    (line 15) //default_random_engine re; // fixed seed

By default, the seed of random_number sequence is automatically changed in each execution (using random_device{}()). Fixed seed will be used if you uncomment 15th line and comment out the 14th line.

### Theory
  This code uses the mathematical theory described in "Filippov, A.V. et al. 2000, *Journal of Colloid and Interface Science* 229, 261–273".

### Contact
nobuhiro.moteki@gmail.com, or, moteki@eps.s.u-tokyo.ac.jp

### Licence
Please see LICENSE.txt
