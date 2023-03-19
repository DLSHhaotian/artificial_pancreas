//
// Created by dlsh on 2021/3/22.
//

#ifndef CASADI_C_TEST_HOVORKAMODEL_H
#define CASADI_C_TEST_HOVORKAMODEL_H
#include <vector>
#include "cmath"
using namespace std;

typedef struct Hovorka_parameters{
    double MwG;
    double AG;
    double tau_D;
    double VG;
    double F01;
    double k12;
    double EGP0;
    double tau_S;
    double VI;
    double ke;
    double ka1;
    double ka2;
    double ka3;
    double SIT;
    double SID;
    double SIE;
    double tau_GI;
    double a;
    double tHR;
    double tin;
    double n;
    double tex;
    double c1;
    double c2;
    double beta;
    double alpha;
    double HRbase;
    double tau_glu;
    double k_glu;
} Hovorka_parameters;

class HovorkaModel {
public:
    HovorkaModel();
    ~HovorkaModel();
    Hovorka_parameters default_parameters();
    vector<double> HovorkaState(double *x, double * u, double d, double HR, double* dw,double* p_sigma);
    double HovarkaOutput(double * x);
    double HovarkaSensor(double * x, double noise);
    int getNstates();
private:
    double MwG;
    double AG;
    double tau_D;
    double VG;
    double F01;
    double k12;
    double EGP0;
    double tau_S;
    double VI;
    double ke;
    double ka1;
    double ka2;
    double ka3;
    double SIT;
    double SID;
    double SIE;
    double tau_GI;
    double a;
    double tHR;
    double tin;
    double n;
    double tex;
    double c1;
    double c2;
    double beta;
    double alpha;
    double HRbase;
    double tau_glu;
    double k_glu;
    int n_states, n_controls, n_disturbance;
    Hovorka_parameters Hovorka_param;
};


#endif //CASADI_C_TEST_HOVORKAMODEL_H
