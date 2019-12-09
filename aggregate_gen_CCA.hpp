
extern std::default_random_engine re;
extern std::uniform_real_distribution<double> urf;

std::vector<Eigen::MatrixXd> aggregate_gen_CCA(const double& a, const int& num_sph_SA, const int& levels, const double& kf, const double& Df, const double& tol);
