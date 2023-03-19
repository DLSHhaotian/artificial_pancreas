//
// Created by dlsh on 2021/3/23.
//
#include "functionUsed.h"
vector<double> functionUsed::NoiseSensorInit(double sigma, double len) {
    double R = sqrt(sigma);
    vector<double> noise(len,0);
    std::normal_distribution<> norm{0, R};
    std::random_device rd;
    std::default_random_engine rng{rd()};
    for(int i=0;i<len;++i){
        noise[i] = norm (rng);
    }
    return noise;
}
vector<double> functionUsed::NoiseSensorUpdate(double noise_cc_prepre, double noise_cc_pre, double noise_v_prepre,
                                               double noise_v_pre, double wcc, double w) {
    vector<double> noise(2,0);
    double noise_cc=1.23*noise_cc_pre-0.3995*noise_cc_prepre+wcc;
    double noise_v=1.013*noise_v_pre-0.2135*noise_v_prepre+w;
    noise[0] = noise_cc;
    noise[1] = noise_v;
    return noise;
}
vector<vector<double>> functionUsed::ScalarStdWienerProcess(double t_final, double N_all, double Ns, double seed) {
    double dt = t_final/N_all;
    double R = sqrt(dt);
    vector<vector<double>> dW(Ns,vector<double>(N_all,0));
    std::normal_distribution<> norm {0, R};
    std::random_device rd;
    std::default_random_engine rng {rd()};
    for(int i=0;i<Ns;++i){
        for(int j=0;j<N_all;++j){
            dW[i][j] = norm (rng);
        }
    }
    return dW;
}