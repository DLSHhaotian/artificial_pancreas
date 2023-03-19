//
// Created by dlsh on 2021/3/22.
//

#include "HovorkaModel.h"
HovorkaModel::HovorkaModel():
        MwG(180.1577),
        AG(0.95273),
        tau_D(50.53245),
        VG(10.602),
        F01(0.7859),
        k12(0.06124),
        EGP0(1.13916),
        tau_S(61.054066),
        VI(9.90742),
        ke(0.174848),
        ka1(0.0077219),
        ka2(0.0911486),
        ka3(0.0273229),
        SIT(59.941792e-4),
        SID(8.40101e-4),
        SIE(227.8809e-4),
        tau_GI(6.7),
        n_states(11),
        n_controls(2),
        n_disturbance(1)
{
    Hovorka_param.MwG = MwG;
    Hovorka_param.AG = AG;
    Hovorka_param.tau_D = tau_D;
    Hovorka_param.VG = VG;
    Hovorka_param.F01 = F01;
    Hovorka_param.k12 = k12;
    Hovorka_param.EGP0 = EGP0;
    Hovorka_param.tau_S = tau_S;
    Hovorka_param.VI = VI;
    Hovorka_param.ke = ke;
    Hovorka_param.ka1 = ka1;
    Hovorka_param.ka2 = ka2;
    Hovorka_param.ka3 = ka3;
    Hovorka_param.SIT = SIT;
    Hovorka_param.SID = SID;
    Hovorka_param.SIE = SIE;
    Hovorka_param.tau_GI = tau_GI;
}
HovorkaModel::~HovorkaModel() {};
Hovorka_parameters HovorkaModel::default_parameters() {
    return Hovorka_param;
}
double HovorkaModel::HovarkaOutput(double *x) {
    return x[4]/VG*18;
}
double HovorkaModel::HovarkaSensor(double *x, double noise) {
    return x[10]*18+noise;
}
vector<double> HovorkaModel::HovorkaState(double *x, double *u, double d, double *dw,
                                          double *p_sigma) {
    vector<double> xdot(n_states,0);
    double uba = u[0], ubo = u[1], dmeal = d;
    double G = x[4]/VG;
    double F01_c = 0, FR = 0;
    if(G>=4.5)
        F01_c = F01;
    else
        F01_c = F01*G/4.5;
    if(G>=9)
        FR = 0.003*(G-9)*VG;
    else
        FR = 0;

    double Dt = 1000*dmeal/MwG, UG = x[1]/tau_D, UI = x[3]/tau_S;
    xdot[0] = AG*Dt-x[0]/tau_D;
    xdot[1] = x[0]/tau_D-x[1]/tau_D;
    xdot[2] = uba+ubo-x[2]/tau_S;
    xdot[3] = x[2]/tau_S-x[3]/tau_S;
    xdot[4] = UG-F01_c-FR-x[7]*x[4]+k12*x[5]+EGP0*(1-x[9]);
    xdot[5] = x[7]*x[4]-(k12+x[8])*x[5];
    xdot[6] = UI/VI-ke*x[6];
    xdot[7] = -ka1*x[7]+SIT*ka1*x[6];
    xdot[8] = -ka2*x[8]+SID*ka2*x[6];
    xdot[9] = -ka3*x[9]+SIE*ka3*x[6];
    xdot[10] = G/tau_GI-x[10]/tau_GI;
    return xdot;

}
int HovorkaModel::getNstates() {
    return n_states;
}