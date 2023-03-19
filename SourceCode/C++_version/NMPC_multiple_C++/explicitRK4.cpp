//
// Created by dlsh on 2021/3/23.
//
#include "explicitRK4.h"
#include <iterator>
using namespace std;
using namespace boost::numeric::odeint;
typedef runge_kutta_dopri5< vector<double> > dopri5_type;
typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;


vector<double> explicitRK4::explicitRK4_Hovorka(double dt, double *xk, double *uk, double dmeal, double HR, double* dw,
                                                double *p_sigma) {
    HovorkaModel Hvok;
    int n_states = Hvok.getNstates();
    double* xk_RK4 = xk;
    double xk_tmp[n_states];
    copy(xk,xk+n_states,xk_tmp);
    int M = 10;
    double DT = dt/M;
    vector<double> RK1, RK2, RK3, RK4;
    for(int i=0;i<M;++i){
        RK1 = Hvok.HovorkaState(xk_RK4,uk,dmeal,HR,dw,p_sigma);
        for(int j=0;j<n_states;++j){
            xk_tmp[j]=xk_RK4[j]+DT/2*RK1[j];
        }
        RK2 = Hvok.HovorkaState(xk_tmp,uk,dmeal,HR,dw,p_sigma);
        for(int j=0;j<n_states;++j){
            xk_tmp[j]=xk_RK4[j]+DT/2*RK2[j];
        }
        RK3 = Hvok.HovorkaState(xk_tmp,uk,dmeal,HR,dw,p_sigma);
        for(int j=0;j<n_states;++j){
            xk_tmp[j]=xk_RK4[j]+DT*RK3[j];
        }
        RK4 = Hvok.HovorkaState(xk_tmp,uk,dmeal,HR,dw,p_sigma);
        for(int j=0;j<n_states;++j){
            xk_RK4[j] = xk_RK4[j]+DT/6*(RK1[j]+2*RK2[j]+2*RK3[j]+RK4[j]);
        }
    }
    vector<double> xk_RK4_vec(xk_RK4,xk_RK4+n_states);
    return xk_RK4_vec;
}
void explicitRK4::write_out(const vector<double> &x, const double t) {
}
vector<double> explicitRK4::Dopri5_Hovorka(double dt, double *xk, double *uk, double dmeal, double HR, double *dw,
                                           double *p_sigma) {
    double err_abs = 1.0e-6;
    double err_rel = 1.0e-3;
    vector<double> xk_vec(xk,xk+16);
    vector<double> u_vec(uk,uk+3), dmeal_vec={dmeal}, HR_vec={HR},p_sigma_vec(p_sigma,p_sigma+16), dw_vec(dw,dw+16);
    vector<vector<double>> params;
    params.push_back(u_vec);
    params.push_back(dmeal_vec);
    params.push_back(HR_vec);
    params.push_back(p_sigma_vec);
    params.push_back(dw_vec);
    struct system mySys(params);
    //controlled_dopri5_type dopri5(make_controlled( err_abs , err_rel , dopri5_type()));
    integrate_adaptive( make_controlled( err_abs , err_rel , dopri5_type()),
                        mySys , xk_vec , 0.0 , dt , 0.5 , write_out);
    return xk_vec;
}
Matrix9d explicitRK4::MVPStateCovar_pred(double *fdot, Matrix9d P_kf, Matrix9d G_mat) {
    Matrix9d fdot_mat,P_pred;
    for(int i=0;i<9;++i){
        for(int j=0;j<9;++j){
            fdot_mat(j,i) = fdot[i*9+j];
        }
    }
    P_pred = fdot_mat*P_kf+P_kf*(fdot_mat.transpose().eval())+G_mat;
    return P_pred;
}

vector<Matrix9d> explicitRK4::explicitRK4_CDEKFpred(double T, double *x, double *u, double d, Matrix9d P_kf, Matrix9d G_mat) {
    MVPmodel MVP;
    int n_statesMVP = 9, M = 4;
    double DT = T/M;
    double xk_RK4[n_statesMVP], xk_tmp[n_statesMVP], fdot[n_statesMVP*n_statesMVP],x_RK1[n_statesMVP],x_RK2[n_statesMVP],x_RK3[n_statesMVP],x_RK4[n_statesMVP];
    copy(x,x+n_statesMVP,xk_RK4);
    copy(x,x+n_statesMVP,xk_tmp);
    Matrix9d P_kf_RK4 = P_kf, P_kf_tmp = P_kf,P_RK1,P_RK2,P_RK3,P_RK4;
    vector<vector<double>> jaco_x;

    for(int k=0;k<M;++k){
        jaco_x = MVP.MVPStateJacobian(xk_RK4,u,d);
        for(int i=0;i<n_statesMVP*n_statesMVP;++i){
            fdot[i] = jaco_x[0][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            x_RK1[i] = jaco_x[1][i];
        }

        P_RK1 = explicitRK4::MVPStateCovar_pred(fdot,P_kf_RK4,G_mat);


        for(int i=0;i<n_statesMVP;++i){
            xk_tmp[i]=xk_RK4[i]+DT/2*x_RK1[i];
        }
        jaco_x = MVP.MVPStateJacobian(xk_tmp,u,d);
        for(int i=0;i<n_statesMVP*n_statesMVP;++i){
            fdot[i] = jaco_x[0][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            x_RK2[i] = jaco_x[1][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            for(int j=0;j<n_statesMVP;++j){
                P_kf_tmp(i,j) = P_kf_RK4(i,j)+DT/2*P_RK1(i,j);
            }
        }
        P_RK2 = explicitRK4::MVPStateCovar_pred(fdot,P_kf_tmp,G_mat);


        for(int i=0;i<n_statesMVP;++i){
            xk_tmp[i]=xk_RK4[i]+DT/2*x_RK2[i];
        }
        jaco_x = MVP.MVPStateJacobian(xk_tmp,u,d);
        for(int i=0;i<n_statesMVP*n_statesMVP;++i){
            fdot[i] = jaco_x[0][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            x_RK3[i] = jaco_x[1][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            for(int j=0;j<n_statesMVP;++j){
                P_kf_tmp(i,j) = P_kf_RK4(i,j)+DT/2*P_RK2(i,j);
            }
        }
        P_RK3 = explicitRK4::MVPStateCovar_pred(fdot,P_kf_tmp,G_mat);


        for(int i=0;i<n_statesMVP;++i){
            xk_tmp[i]=xk_RK4[i]+DT*x_RK3[i];
        }
        jaco_x = MVP.MVPStateJacobian(xk_tmp,u,d);
        for(int i=0;i<n_statesMVP*n_statesMVP;++i){
            fdot[i] = jaco_x[0][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            x_RK4[i] = jaco_x[1][i];
        }
        for(int i=0;i<n_statesMVP;++i){
            for(int j=0;j<n_statesMVP;++j){
                P_kf_tmp(i,j) = P_kf_RK4(i,j)+DT*P_RK3(i,j);
            }
        }
        P_RK4 = explicitRK4::MVPStateCovar_pred(fdot,P_kf_tmp,G_mat);


        for(int i=0;i<n_statesMVP;++i){
            xk_RK4[i]=xk_RK4[i]+DT/6*(x_RK1[i]+2*x_RK2[i]+2*x_RK3[i]+x_RK4[i]);
        }

        for(int i=0;i<n_statesMVP;++i){
            for(int j=0;j<n_statesMVP;++j){
                P_kf_RK4(i,j) = P_kf_RK4(i,j)+DT/6*(P_RK1(i,j)+2*P_RK2(i,j)+2*P_RK3(i,j)+P_RK4(i,j));
            }
        }

    }
    vector<Matrix9d> xk_Pk;
    Matrix9d xk_mat;
    for(int i=0;i<n_statesMVP;++i){
        xk_mat(i,0) = xk_RK4[i];
    }
    xk_Pk.push_back(xk_mat);
    xk_Pk.push_back(P_kf_RK4);
    return xk_Pk;
}