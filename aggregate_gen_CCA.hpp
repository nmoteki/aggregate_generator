
extern std::default_random_engine re;
extern std::uniform_real_distribution<float> urf;

Eigen::MatrixXf aggregate_gen_CCA(const float& a, const int& num_sph_SA, const int& levels, const float& kf, const float& Df, const float& tol);
