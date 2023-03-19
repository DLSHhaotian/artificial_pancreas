//
// Created by dlsh on 2021/3/18.
//

#ifndef CASADI_C_TEST_MVPMODEL_H
#define CASADI_C_TEST_MVPMODEL_H
#include <iostream>
#include <vector>
using namespace std;

typedef struct MVP_parameters{
    double k1;
    double CI;
    double SI;
    double EGP;
    double km;
    double VG;
    double tau_GI;
    double sigma4;
    double tau_glu;
    double k_glu;
} MVP_parameters;

class MVPmodel {
public:
    MVPmodel();
    ~MVPmodel();
    MVP_parameters default_parameters();
    vector<vector<double>> MVPOutputJacobian(double* x);
    vector<vector<double>> MVPStateJacobian(double* x, double* u, double d);
    //[hdot,z] = MVPOutputJacobian(x)
    //[fdot,xdot] = MVPStateJacobian(x,u,d,p)
private:
    double k1;
    double CI;
    double SI;
    double EGP;
    double km;
    double VG;
    double tau_GI;
    double sigma4;
    double tau_glu;
    double k_glu;
    int n_states, n_controls, n_disturbance;
    MVP_parameters MVP_param;
};


#endif //CASADI_C_TEST_MVPMODEL_H
