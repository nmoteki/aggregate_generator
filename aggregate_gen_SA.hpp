
extern std::default_random_engine re;
extern std::uniform_real_distribution<float> urf;

Eigen::MatrixXf aggregate_gen_SA(const float& a, const int& num_sph, const float& kf, const float& Df, const float& tol);
