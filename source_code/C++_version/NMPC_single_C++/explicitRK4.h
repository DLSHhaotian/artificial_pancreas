//
// Created by dlsh on 2021/3/23.
//

#ifndef CASADI_C_TEST_EXPLICITRK4_H
#define CASADI_C_TEST_EXPLICITRK4_H
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

#include "HovorkaModel.h"
#include "MVPmodel.h"
using namespace std;
using namespace Eigen;
typedef Eigen::Matrix<double,7,7> Matrix7d;

struct system
{
    vector<vector<double>> sys_param;
    system( vector<vector<double>> param ) : sys_param( param ) {}

    void operator()( const vector<double> &x , vector<double> &dxdt , const double t ) const
    {
        const double MwG = 180.1577, AG = 0.95273, tau_D = 50.53245, VG = 10.602, F01 = 0.7859, k12 = 0.06124, EGP0 = 1.13916, tau_S = 61.054066,VI = 9.90742,ke = 0.174848, ka1 = 0.0077219, ka2 = 0.0911486, ka3= 0.0273229, SIT = 59.941792e-4,SID = 8.40101e-4,SIE = 227.8809e-4,tau_GI = 6.7;
        vector<double> u = sys_param[0], d = sys_param[1], p_sigma = sys_param[2], dw = sys_param[3];
        double uba = u[0], ubo = u[1], dmeal = d[0];
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
        dxdt[0] = AG*Dt-x[0]/tau_D;
        dxdt[1] = x[0]/tau_D-x[1]/tau_D;
        dxdt[2] = uba+ubo-x[2]/tau_S;
        dxdt[3] = x[2]/tau_S-x[3]/tau_S;
        dxdt[4] = UG-F01_c-FR-x[7]*x[4]+k12*x[5]+EGP0*(1-x[9]);
        dxdt[5] = x[7]*x[4]-(k12+x[8])*x[5];
        dxdt[6] = UI/VI-ke*x[6];
        dxdt[7] = -ka1*x[7]+SIT*ka1*x[6];
        dxdt[8] = -ka2*x[8]+SID*ka2*x[6];
        dxdt[9] = -ka3*x[9]+SIE*ka3*x[6];
        dxdt[10] = G/tau_GI-x[10]/tau_GI;
    }
};
namespace explicitRK4{
    vector<double> explicitRK4_Hovorka(double dt,double* xk,double* uk,double dmeal,double* dw,double* p_sigma);
    vector<double> Dopri5_Hovorka(double dt,double* xk,double* uk,double dmeal,double* dw,double* p_sigma);
    void write_out( const vector<double> &x , const double t );
    Matrix7d MVPStateCovar_pred(double* fdot,Matrix7d P_kf,Matrix7d G_mat);
    vector<Matrix7d> explicitRK4_CDEKFpred(double T, double* x, double* u, double d, Matrix7d P_kf, Matrix7d G_mat);

}
#endif //CASADI_C_TEST_EXPLICITRK4_H
