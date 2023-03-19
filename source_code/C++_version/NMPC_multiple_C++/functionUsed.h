//
// Created by dlsh on 2021/3/23.
//

#ifndef CASADI_C_TEST_FUNCTIONUSED_H
#define CASADI_C_TEST_FUNCTIONUSED_H
#include <vector>
#include <random>
using namespace std;
namespace functionUsed{
    vector<double> NoiseSensorInit(double sigma, double len);
    vector<double> NoiseSensorUpdate(double noise_cc_prepre,double noise_cc_pre,double noise_v_prepre,double noise_v_pre,double wcc,double w);
    vector<vector<double>> ScalarStdWienerProcess(double t_final, double N_all,double Ns,double seed);
}
#endif //CASADI_C_TEST_FUNCTIONUSED_H
