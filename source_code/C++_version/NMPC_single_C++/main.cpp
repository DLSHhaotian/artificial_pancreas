#include <iostream>
#include <fstream>
#include <ctime>
#include <casadi/casadi.hpp>
#include <vector>
#include <string>
#include <map>
#include "MVPmodel.h"
#include "HovorkaModel.h"
#include "functionUsed.h"
#include "explicitRK4.h"
#include "CasADiSolver.h"
#include <float.h>
#include <time.h>
#include <Eigen/Dense>
#include <random>
#include <fstream>

using namespace Eigen;
using namespace casadi;
using namespace std;

typedef Eigen::Matrix<double,7,7> Matrix7d;
typedef Eigen::Matrix<double,1,7> Matrix1_7d;
typedef Eigen::Matrix<double,7,1> Matrix7_1d;

int main() {
    const int t_final = 24*60, T = 5, N = 72;
    const int N_all = t_final/T+1;
    const int n_statesHovk = 11,n_statesMVP =7,n_controls = 2,n_disturbance = 1;
    double z_ref = 108;
    double x0_Hovk[n_statesHovk] = {0.0,0.0,664.58396,664.58396,63.61117887,36.02065,6.28366995,0.0376654,0.0052789,0.14319,6.0};
    double x0_MVP[n_statesMVP] = {10.043388,10.043388, 0.01999, 108.0000, 0, 0, 108.0000};
    double t_arr[N_all], x_arr[n_statesHovk*N_all],y_arr[N_all],z_arr[N_all],u_arr[(N_all-1)*n_controls],d_arr[(N_all-1)*n_disturbance],t_cost[N_all-1];
    double x_bar[n_statesMVP*N_all],z_bar[N_all];

    for(int i=0;i<N_all;++i){
        t_arr[i] = i*T;
    }
    for(int i=0;i<N_all-1;++i){
        d_arr[i] = 0;
    }
    d_arr[6*60/T-1] = 50/T;
    //d_arr[0] = 50/T;
    d_arr[12*60/T-1] = 60/T;
    d_arr[18*60/T-1] = 70/T;

    for(int i=0;i<n_statesHovk;++i){
        x_arr[i]=x0_Hovk[i];
    }

    for(int i=0;i<n_statesMVP;++i){
        x_bar[i]=x0_MVP[i];
    }

    //process noise initialization

    double Ns = n_statesHovk,seed =100;
    vector<vector<double>> dW =functionUsed::ScalarStdWienerProcess(t_final,N_all,Ns,seed);
    double dW_k[n_statesHovk], p_sigma[n_statesHovk];
    for(int i = 0;i<n_statesHovk;++i){
        p_sigma[i] = x0_Hovk[i]*0.01;
    }
    //sensor noise initialization
    double R_wcc = 11.3, R_w =14.45;
    double noise_cc_prepre=0, noise_cc_pre=0, noise_v_prepre=0, noise_v_pre=0,noise_cc_k = 0,noise_v_k = 0;
    double wcc_arr[N_all],w_arr[N_all];
    vector<double> wcc_arr_vec = functionUsed::NoiseSensorInit(R_wcc,N_all);
    vector<double> w_arr_vec = functionUsed::NoiseSensorInit(R_w,N_all);
    for(int i=0;i<N_all;++i){
        wcc_arr[i] = wcc_arr_vec[i];       w_arr[i] = w_arr_vec[i];
    }


    //CDEKF initialization
    double sigma4 = 3.05, e_kf = 0, R_kf = 0.02*0.02;
    Matrix1_7d hdot;
    Matrix7_1d K_kf;
    vector<vector<double>> jaco_z;
    vector<Matrix7d> pred_x;
    Matrix7d P_kf;
    for(int i=0;i<n_statesMVP;++i){
        for(int j=0;j<n_statesMVP;++j){
            P_kf(i,j) = 0;
        }
    }
    P_kf(3,3) = sigma4*sigma4;
    Matrix7d G_mat = P_kf;

    //initial guess of decision variables
    double OPT_varX0[n_statesMVP*(N+1)+n_controls*N];
    for(int i=0;i<N+1;++i){
        for(int j=0;j<n_statesMVP;++j){
            OPT_varX0[n_statesMVP*i+j] = x0_MVP[j];
        }
    }
    for(int i=n_statesMVP*(N+1);i<n_statesMVP*(N+1)+n_controls*N;++i){
        OPT_varX0[i] = 0;
    }
    vector<double> OPT_varX0_debug(OPT_varX0,OPT_varX0+n_statesMVP*(N+1)+n_controls*N);

    HovorkaModel Hvok;
    MVPmodel MVP;
    CasADiSolver nlpsolver;
    map<string, vector<double>> args;
    Function solver;
    vector<double> sol,xsol_debug,usol_debug;
    double xsol[n_statesMVP*(N+1)],usol[n_controls*N];

    vector<double> xk_vec(n_statesHovk,0);
    double yk = 0, zk = 0, zk_bar = 0, xk[n_statesHovk],xk_bar[n_statesMVP] ,uk[2] = {0,0},uk_d[3] = {0,0,0};


    clock_t c_start, c_end;
    //cout<<"-DBL_MAX:"<<-DBL_MAX<<endl;
    //cout<<"DBL_MIN:"<<DBL_MIN<<endl;
    int k = 0;
    //for(;k<1;++k){
    for(;k<N_all-1;++k){
        cout<<"The iterator:"<<k<<endl;
        vector<double> noise_k = functionUsed::NoiseSensorUpdate(noise_cc_prepre,noise_cc_pre,noise_v_prepre,noise_v_pre,wcc_arr[k],w_arr[k]);
        noise_cc_k = noise_k[0];
        noise_v_k = noise_k[1];

        for(int i = 0;i<n_statesHovk;++i){
            xk[i] = x_arr[k*n_statesHovk+i];
        }
        yk = Hvok.HovarkaSensor(xk,0);
        zk = Hvok.HovarkaOutput(xk);
        y_arr[k] = yk;
        z_arr[k] = zk;
        cout<<"The output end:"<<k<<endl;

        c_start = clock();
        //CDEKF predict
        for(int i = 0;i<n_statesMVP;++i){
            xk_bar[i] = x_bar[k*n_statesMVP+i];
        }
        uk_d[0] = uk[0], uk_d[1] = uk[1], uk_d[2] = d_arr[k]*1000;
        for(int i=0;i<n_statesHovk;++i){
            //dW_k[i] = dW[i][k];
            dW_k[i] = 0;
        }
        R_kf = 0.02*0.02;
        pred_x = explicitRK4::explicitRK4_CDEKFpred(T,xk_bar,uk_d,dW_k[0],P_kf,G_mat);
        Matrix7d xk_bar_mat = pred_x[0];
        for(int i=0;i<n_statesMVP;++i){
            xk_bar[i] = xk_bar_mat(i,0);
        }
        P_kf = pred_x[1];
        cout<<"The CDEKF predict end:"<<k<<endl;
        //update
        jaco_z = MVP.MVPOutputJacobian(xk_bar);
        for(int i=0;i<n_statesMVP;++i){
            hdot(0,i) = jaco_z[0][i];
        }
        zk_bar = jaco_z[1][0];
        e_kf = y_arr[k]-zk_bar;
        R_kf = hdot*P_kf*(hdot.transpose().eval())+R_kf;
        K_kf = P_kf*(hdot.transpose().eval())/R_kf;
        Matrix7_1d xk_bar_plus = K_kf*e_kf;
        for(int i=0;i<n_statesMVP;++i){
            xk_bar[i] = xk_bar[i]+xk_bar_plus(i,0);
        }
        P_kf = P_kf-K_kf*R_kf*(K_kf.transpose().eval());
        for(int i=0;i<n_statesMVP;++i){
            x_bar[(k+1)*n_statesMVP+i]=xk_bar[i];
        }
        z_bar[k] = zk_bar;
        cout<<"The CDEKF update end:"<<k<<endl;
        //cout<<"xk_bar:"<<xk_bar[0]<<","<<xk_bar[1]<<","<<xk_bar[2]<<","<<xk_bar[3]<<","<<xk_bar[4]<<","<<xk_bar[5]<<","<<xk_bar[6]<<endl;

        nlpsolver.NLPsolverConfig(OPT_varX0,xk_bar,d_arr[k]*1000,z_ref);
        args = nlpsolver.getArgs();
        solver = nlpsolver.getSolver();
        cout<<"The nlpsolver end:"<<k<<endl;
        solver({{"lbx", args["lbx"]}, {"ubx", args["ubx"]}, {"x0", args["x0"]}, {"lbg", args["lbg"]}, {"ubg", args["ubg"]},{"p",args["p"]}}, {{"x", &sol}});
        c_end   = clock();
        //cout<<"sol"<<sol<<endl;
        printf("The pause %d used %f ms by clock()\n",k,difftime(c_end,c_start)/1000);
        t_cost[k] = difftime(c_end,c_start)/1000000;
        for(int i=0;i<n_statesMVP*(N+1);++i){
            xsol[i] = sol[i];
        }
        for(int i=0;i<n_controls*N;++i){
            usol[i] = sol[n_statesMVP*(N+1)+i];
        }
        uk[0] = usol[0], uk[1] = usol[1];
        for(int i=0;i<n_controls;++i){
            u_arr[k*n_controls+i]=uk[i];
        }

        vector<double> xsol_debug(xsol,xsol+n_statesMVP*(N+1)),usol_debug(usol,usol+n_controls*N);
        if(d_arr[k]>0){
            cout<<"The xsol:"<<xsol_debug<<endl;
            cout<<"The usol:"<<usol_debug<<endl;
            cout<<"The uk:"<<uk[0]<<", "<<uk[1]<<endl;
        }
        sol.clear();
        for(int i=0;i<n_statesMVP*N;++i){
            OPT_varX0[i] = xsol[i+n_statesMVP];
        }
        for(int i=0;i<n_statesMVP;++i){
            OPT_varX0[n_statesMVP*N+i] = xsol[n_statesMVP*N+i];
        }
        for(int i=0;i<n_controls*(N-1);++i){
            OPT_varX0[n_statesMVP*(N+1)+i] = usol[i+n_controls];
        }
        for(int i=0;i<n_controls;++i){
            OPT_varX0[n_statesMVP*(N+1)+n_controls*(N-1)+i] = usol[n_controls*(N-1)+i];
        }
        cout<<"The init guess end:"<<k<<endl;

        //cout<<"p_sigma*dW_k:"<<p_sigma[0]*dW_k[0]<<" "<<p_sigma[1]*dW_k[1]<<" "<<p_sigma[2]*dW_k[2]<<" "<<p_sigma[3]*dW_k[3]<<" "<<p_sigma[4]*dW_k[4]<<" "<<p_sigma[5]*dW_k[5]<<" "<<p_sigma[6]*dW_k[6]<<" "<<p_sigma[7]*dW_k[7]<<" "<<p_sigma[8]*dW_k[8]<<" "<<p_sigma[9]*dW_k[9]<<" "<<p_sigma[10]*dW_k[10]<<endl;
        //xk_vec = explicitRK4::explicitRK4_Hovorka(T,xk,uk,d_arr[k],dW_k,p_sigma);
        xk_vec = explicitRK4::Dopri5_Hovorka(T,xk,uk,d_arr[k],dW_k,p_sigma);
        cout<<"The ode end:"<<k<<endl;
        for(int i = 0;i<n_statesHovk;++i){
            x_arr[(k+1)*n_statesHovk+i] = xk_vec[i];
        }
    }
    k = N_all-1;
    yk = Hvok.HovarkaSensor(xk,noise_cc_k+noise_v_k);
    zk = Hvok.HovarkaOutput(xk);
    y_arr[k] = yk;
    z_arr[k] = zk;
    z_bar[k] = zk_bar;


    ofstream outfile_t_arr,outfile_x_arr,outfile_u_arr,outfile_d_arr,outfile_y_arr,outfile_z_arr,outfile_t_cost,outfile_x_bar,outfile_z_bar;
    outfile_t_arr.open("t_arr.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_x_arr.open("x_arr.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_u_arr.open("u_arr.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_d_arr.open("d_arr.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_y_arr.open("y_arr.txt", ios::binary | ios::trunc| ios::in | ios::out);
    outfile_z_arr.open("z_arr.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_t_cost.open("t_cost.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_x_bar.open("x_bar.txt", ios::binary | ios::trunc | ios::in | ios::out);
    outfile_z_bar.open("z_bar.txt", ios::binary | ios::trunc | ios::in | ios::out);
    for(int i=0;i<N_all;++i){
        outfile_t_arr<<t_arr[i]<<"\n";
        outfile_x_arr<<x_arr[i*n_statesHovk+0]<<" "<<x_arr[i*n_statesHovk+1] <<" "<< x_arr[i*n_statesHovk+2]<<" "<<x_arr[i*n_statesHovk+3] <<" "<<x_arr[i*n_statesHovk+4]<<" "<<x_arr[i*n_statesHovk+5] <<" "<<x_arr[i*n_statesHovk+6]<<" "<<x_arr[i*n_statesHovk+7] <<" "<<x_arr[i*n_statesHovk+8]<<" "<<x_arr[i*n_statesHovk+9] <<" "<<x_arr[i*n_statesHovk+10]<<"\n";
        outfile_y_arr<<y_arr[i]<<"\n";
        outfile_z_arr<<z_arr[i]<<"\n";
        outfile_z_bar<<z_bar[i]<<"\n";
        outfile_x_bar<<x_bar[i*n_statesMVP+0]<<" "<<x_bar[i*n_statesMVP+1] <<" "<< x_bar[i*n_statesMVP+2]<<" "<<x_bar[i*n_statesMVP+3] <<" "<<x_bar[i*n_statesMVP+4]<<" "<<x_bar[i*n_statesMVP+5] <<" "<<x_bar[i*n_statesMVP+6]<<"\n";
        if(i!=(N_all-1)){
            outfile_u_arr<<u_arr[i*n_controls+0]<<" "<<u_arr[i*n_controls+1]<<"\n";
            outfile_d_arr<<d_arr[i]<<"\n";
            outfile_t_cost<<t_cost[i]<<"\n";
        }
    }

    outfile_t_arr.close();
    outfile_x_arr.close();
    outfile_u_arr.close();
    outfile_d_arr.close();
    outfile_y_arr.close();
    outfile_z_arr.close();
    outfile_t_cost.close();
    outfile_x_bar.close();
    outfile_z_bar.close();

    return 0;
}
/*
MVPmodel MVP;
MVP_parameters MVP_param = MVP.default_parameters();
//p = [k1; CI; SI; EGP; km; VG; tau_GI; sigma4];
SX x1 = SX::sym("x1"); SX x2 = SX::sym("x2"); SX x3 = SX::sym("x3");
SX x4 = SX::sym("x4"); SX x5 = SX::sym("x5"); SX x6 = SX::sym("x6"); SX x7 = SX::sym("x7");
vector<SX> states = {x1,x2,x3,x4,x5,x6,x7};
int n_states = states.size();
SX states_arr = vertcat(states);


SX u1 = SX::sym("u1"); SX u2 = SX::sym("u2");
vector<SX> controls = {u1,u2};
int n_controls = controls.size();
SX controls_arr = vertcat(controls);

SX d = SX::sym("d");
vector<SX> disturbance = {d};
SX disturbance_arr = vertcat(disturbance);


vector<SX> xdot(7);
xdot[0] = MVP_param.k1*((u1+u2)/MVP_param.CI-x1);
xdot[1] = MVP_param.k1*(x1-x2);
xdot[2] = -MVP_param.k1*x3+MVP_param.k1*MVP_param.SI*x2;
xdot[3] = -x3*x4+MVP_param.EGP+(x6*MVP_param.km/MVP_param.VG);
xdot[4] = d-x5*MVP_param.km;
xdot[5] = (x5-x6)*MVP_param.km;
xdot[6] = (x4-x7)/MVP_param.tau_GI;

vector<SX> XUD(1,states_arr);
XUD.push_back(controls_arr);
XUD.push_back(disturbance_arr);
Function f = Function("f", XUD,xdot);

cout << "f: " << f<< endl;
//ERK4 of numerical method
const double T = 5;//[min]
const int M = 4;
double DT = T/M;
SX states_arr_RK4 = states_arr;
SX states_arr_tmp = states_arr;
vector<SX> XUD_RK4;
vector<SX> RK1(1),RK2(1),RK3(1),RK4(1);
for(int i=0;i<M;++i){

XUD_RK4.push_back(states_arr_RK4);
XUD_RK4.push_back(controls_arr);
XUD_RK4.push_back(disturbance_arr);
RK1 = f(XUD_RK4);
XUD_RK4.clear();
cout << "RK1: " << RK1 << endl;

states_arr_tmp = states_arr_RK4+DT/2*vertcat(RK1);
XUD_RK4.push_back(states_arr_tmp);
XUD_RK4.push_back(controls_arr);
XUD_RK4.push_back(disturbance_arr);
RK2 = f(XUD_RK4);
XUD_RK4.clear();

states_arr_tmp = states_arr_RK4+DT/2*vertcat(RK2);
XUD_RK4.push_back(states_arr_tmp);
XUD_RK4.push_back(controls_arr);
XUD_RK4.push_back(disturbance_arr);
RK3 = f(XUD_RK4);
XUD_RK4.clear();

states_arr_tmp = states_arr_RK4+DT*vertcat(RK3);
XUD_RK4.push_back(states_arr_tmp);
XUD_RK4.push_back(controls_arr);
XUD_RK4.push_back(disturbance_arr);
RK4 = f(XUD_RK4);
XUD_RK4.clear();

states_arr_RK4 = states_arr_RK4 + DT/6*(vertcat(RK1)+2*vertcat(RK2)+2*vertcat(RK3)+vertcat(RK4));
}
Function fun_RK4("F", XUD,vector<SX>{states_arr_RK4});
cout << "fun_RK4: " << fun_RK4 << endl;
//Self-defined functions
SX val_max = log(exp(x1)+exp(x2));
SX val_min = -log(exp(-x1)+exp(-x2));
SX val_abs = sqrt(x1*x1+1);
Function myMax("myMax",{x1,x2},{val_max});
Function myMin("myMin",{x1,x2},{val_min});
Function myAbs("myAbs",{x1},{val_abs});

const int N =72;

MX X = MX::sym("X",n_states*(N+1));
MX U = MX::sym("U",n_controls*N);
MX P = MX::sym("P",n_states+1+1);
MX OPT_variables = vertcat(X,U);

vector<MX> X_vec(1,X(Slice(0,n_states),Slice()));
vector<MX> U_vec(1,U(Slice(0,n_controls),Slice()));
for(int i=1;i<N+1;++i){
X_vec.push_back(X(Slice(i*n_states,(i+1)*n_states),Slice()));
if(i!=N){
U_vec.push_back(U(Slice(i*n_controls,(i+1)*n_controls),Slice()));
}
}
const double Qz = 0.01, Qz_max = 1, Qz_min = 300000, Qu_bo = 100, Qu_delta_ba = 150;
const double z_max = 180, z_min = 70;
MX obj = 0;
for(int i=0;i<N;++i){
MX st = X_vec[i];
MX zk = st(3), con = U_vec[i];
MX max_k = vertcat(myMax({zk - z_max, 0}));
MX min_k = vertcat(myMin({zk - z_min, 0}));
MX abs_k = vertcat(myAbs(con(1)));
obj = obj + (zk-P(n_states))*Qz*(zk-P(n_states)) + max_k*Qz_max*max_k + min_k*Qz_min*min_k + Qu_bo*abs_k;
if(i!=N-1){
MX con_next = U_vec[i+1];
obj = obj + (con_next(0)-con(0))*Qu_delta_ba*(con_next(0)-con(0));
}
}

vector<MX> g;
MX st = X_vec[0];
MX zk = st(3);
g.push_back(st-P(Slice(0,n_states),Slice()));
for(int i=0;i<N;++i){
st = X_vec[i];
zk = st(3);
MX st_next = X_vec[i+1], con = U_vec[i], st_next_RK4;
if(i==0){
vector<MX> XUD_tmp(1,st);
XUD_tmp.push_back(con);
XUD_tmp.push_back(P(n_states+1));
st_next_RK4 = vertcat(fun_RK4(XUD_tmp));
}
else{
vector<MX> XUD_tmp(1,st);
XUD_tmp.push_back(con);
XUD_tmp.push_back(0);
st_next_RK4 = vertcat(fun_RK4(XUD_tmp));
}

g.push_back(st_next-st_next_RK4);
}


MXDict nlp_prob = {{"x", OPT_variables}, {"f", obj}, {"g", vertcat(g)}, {"p",P}};
Dict opts;
opts["ipopt.tol"] = 1e-8;
opts["ipopt.max_iter"] = 500;
opts["ipopt.print_level"] = 5;
Function solver = nlpsol("nlpsol", "ipopt", nlp_prob, opts);

map<string, vector<double>> args;
//map<string, DM> args,res;
vector<double> x0 = {12.5805,12.5805, 0.0116, 108.0000, 0, 0, 108.0000};
vector<double> OPT_varX0,OPT_varU0(n_controls*N,0);
for(int i=0;i<N+1;++i){
OPT_varX0.insert(OPT_varX0.end(),x0.begin(),x0.end());
}
OPT_varX0.insert(OPT_varX0.end(),OPT_varU0.begin(),OPT_varU0.end());


const double u_min = 0, u_max = DBL_MAX, x_min = -DBL_MAX, x_max = DBL_MAX;
vector<double> OPT_varX_min(n_states*(N+1),x_min);
vector<double> OPT_varX_max(n_states*(N+1),x_max);
vector<double> OPT_varU_min(n_controls*N,u_min);
vector<double> OPT_varU_max(n_controls*N,u_max);
OPT_varX_min.insert(OPT_varX_min.end(),OPT_varU_min.begin(),OPT_varU_min.end());
OPT_varX_max.insert(OPT_varX_max.end(),OPT_varU_max.begin(),OPT_varU_max.end());

vector<double> x_bar = x0;
double z_ref = 108, d_meal = 50;
vector<double> OPT_p = x_bar;
OPT_p.push_back(z_ref);
OPT_p.push_back(d_meal);

args["lbx"] = OPT_varX_min;
args["ubx"] = OPT_varX_max;
args["lbg"] = vector<double>(n_states*(N+1),0);
args["ubg"] = vector<double>(n_states*(N+1),0);
args["x0"] = OPT_varX0;
args["p"] = OPT_p;


clock_t c_start, c_end;
c_start = clock();



cout << "OPT_varX_min: " << OPT_varX_min << endl;
cout << "obj: " << obj << endl;
cout << "g: " << g << endl;
//res = solver(args);


c_end   = clock();
printf("The pause used %f ms by clock()\n",difftime(c_end,c_start)/1000);
vector<double> sol,xsol,usol,uk;
solver({{"lbx", args["lbx"]}, {"ubx", args["ubx"]}, {"x0", args["x0"]}, {"lbg", args["lbg"]}, {"ubg", args["ubg"]},{"p",args["p"]}}, {{"x", &sol}});
//solver({{"lbx", OPT_varX_min}, {"ubx", OPT_varX_max}, {"x0", OPT_varX0}, {"lbg", vector<double>(N+1,0)}, {"ubg", vector<double>(N+1,0)},{"p",OPT_p}}, {{"x", &Usol}});
//cout << "Usol: " << Usol << endl;
xsol.insert(xsol.begin(),sol.begin(),sol.begin()+n_states*(N+1));
usol.insert(usol.begin(),sol.begin()+n_states*(N+1),sol.end());
uk.insert(uk.begin(),usol.begin(),usol.begin()+2);
OPT_varX0.clear();
OPT_varX0.insert(OPT_varX0.begin(),xsol.begin()+n_states,xsol.end());
OPT_varX0.insert(OPT_varX0.end(),xsol.begin()+n_states*N,xsol.end());
OPT_varX0.insert(OPT_varX0.end(),usol.begin()+n_controls,usol.end());
OPT_varX0.insert(OPT_varX0.end(),usol.begin()+n_controls*(N-1),usol.end());

cout << "xsol: " << xsol << endl;
cout << "usol: " << usol << endl;
cout << "OPT_varX0: " << OPT_varX0 << endl;
cout << "OPT_varX0 size: " << OPT_varX0.size() << endl;
*/



//std::cout << sol.value(x) << ":" << sol.value(y) << std::endl;