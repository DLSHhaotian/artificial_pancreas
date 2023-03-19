//
// Created by dlsh on 2021/3/23.
//

#ifndef CASADI_C_TEST_CASADISOLVER_H
#define CASADI_C_TEST_CASADISOLVER_H
#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include "MVPmodel.h"
#include <float.h>
using namespace std;
using namespace casadi;
class CasADiSolver {
public:
    CasADiSolver();
    ~CasADiSolver();
    void NLPsolverConfig(double* x0,double* xbar,double dmeal,double z_ref);
    map<string, vector<double>> getArgs();
    Function getSolver();
private:
    MVPmodel MVP;
    const int T = 5, M = 5, N = 72, n_states = 9, n_controls = 3, n_disturbance = 1;
    const double Qz = 0.01, Qz_max = 1, Qz_min = 300000, Qu_bo = 30, Qu_delta_ba = 150, Qu_glu = 300;
    const double z_max = 180, z_min = 80, uba_ref = 15.0651;
    const double u_min = 0, u_max = DBL_MAX, x_min = -DBL_MAX, x_max = DBL_MAX;
    SX x1,x2,x3,x4,x5,x6,x7,x8,x9;
    SX u1,u2,u3;
    SX d;
    SX states_arr,controls_arr,disturbance_arr;
    Function f,fun_RK4,myMax,myMin,myAbs;
    MX X,U,P,OPT_variables,obj,g;
    MXDict nlp_prob;
    Dict opts;
    Function solver;
    map<string, vector<double>> args;
};
#endif //CASADI_C_TEST_CASADISOLVER_H
