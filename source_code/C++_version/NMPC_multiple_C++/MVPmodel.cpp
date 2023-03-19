//
// Created by dlsh on 2021/3/18.
//

#include "MVPmodel.h"
//tau_glu=19;%[min] glucagon
//        k_glu=0.075;
MVPmodel::MVPmodel():
    k1(0.011365),
    CI(1.5),
    SI(0.0019905),
    EGP(2.15909),
    km(0.01347),
    VG(119.14286),
    tau_GI(6.7),
    sigma4(0.344495),
    tau_glu(19),
    k_glu(0.075),
    n_states(9),
    n_controls(3),
    n_disturbance(1)
{
    MVP_param.k1 = k1;
    MVP_param.CI = CI;
    MVP_param.SI = SI;
    MVP_param.EGP = EGP;
    MVP_param.km = km;
    MVP_param.VG = VG;
    MVP_param.tau_GI = tau_GI;
    MVP_param.sigma4 = sigma4;
    MVP_param.tau_glu = tau_glu;
    MVP_param.k_glu = k_glu;
}
MVPmodel::~MVPmodel() {}
MVP_parameters MVPmodel::default_parameters() {
    return MVP_param;
}
vector<vector<double>> MVPmodel::MVPOutputJacobian(double* x) {
    vector<vector<double>> hdot_z;
    vector<double> hdot(n_states,0);
    vector<double> z(1,x[3]);
    hdot[3] = 1;
    hdot_z.push_back(hdot);
    hdot_z.push_back(z);
    return hdot_z;
}

vector<vector<double>> MVPmodel::MVPStateJacobian(double* x, double* u, double d) {
    vector<vector<double>> fdot_xdot;
    vector<double> fdot(n_states*n_states,0);
    vector<double> xdot(n_states,0);
    //MVP p=[k1; CI; SI; EGP; km; VG; tau_GI; rou4]
    fdot[0] = -k1;
    fdot[1] = k1;
    fdot[1*n_states+1] = -k1;
    fdot[1*n_states+2] = SI*k1;
    fdot[2*n_states+2] = -k1;
    fdot[2*n_states+3] = -x[3];
    fdot[3*n_states+3] = -x[2];
    fdot[5*n_states+3] = km/VG;
    fdot[8*n_states+3] = k_glu;
    fdot[4*n_states+4] = -1/tau_GI;
    fdot[4*n_states+5] = km;
    fdot[5*n_states+5] = -km;
    fdot[3*n_states+6] = 1/tau_GI;
    fdot[6*n_states+6] = -1/tau_GI;
    fdot[7*n_states+7] = -1/tau_glu;
    fdot[7*n_states+8] = 1/tau_glu;
    fdot[8*n_states+8] = -1/tau_glu;
    double uba = u[0], ubo = u[1], uglu = u[2], dmeal = u[3];
    double dw = d;
    xdot[0] = k1*((uba+ubo)/CI - x[0]);
    xdot[1] = k1*(x[0]-x[1]);
    xdot[2] = -k1*x[2]+k1*SI*x[1];
    xdot[3] = -x[2]*x[3]+EGP+(x[5]*km/VG)+k_glu*x[8]+sigma4*dw;
    xdot[4] = dmeal-x[4]*km;
    xdot[5] = (x[4]-x[5])*km;
    xdot[6] = (x[3]-x[6])/tau_GI;
    xdot[7] = uglu-x[7]/tau_glu;
    xdot[8] = (x[7]-x[8])/tau_glu;
    //cout<<"xdot:"<< xdot[4] <<endl;

    fdot_xdot.push_back(fdot);
    fdot_xdot.push_back(xdot);
    return fdot_xdot;
}