//
// Created by dlsh on 2021/3/23.
//
#include "CasADiSolver.h"
CasADiSolver::CasADiSolver(){
    MVP_parameters MVP_param = MVP.default_parameters();
    x1 = SX::sym("x1"); x2 = SX::sym("x2"); x3 = SX::sym("x3");
    x4 = SX::sym("x4"); x5 = SX::sym("x5"); x6 = SX::sym("x6");
    x7 = SX::sym("x7"); x8 = SX::sym("x8"); x9 = SX::sym("x9");
    vector<SX> states = {x1,x2,x3,x4,x5,x6,x7,x8,x9};
    int n_states = states.size();
    states_arr = vertcat(states);

    u1 = SX::sym("u1"); u2 = SX::sym("u2"); u3 = SX::sym("u3");
    vector<SX> controls = {u1,u2,u3};
    int n_controls = controls.size();
    controls_arr = vertcat(controls);

    d = SX::sym("d");
    vector<SX> disturbance = {d};
    disturbance_arr = vertcat(disturbance);

    vector<SX> xdot(9);
    xdot[0] = MVP_param.k1*((u1+u2)/MVP_param.CI-x1);
    xdot[1] = MVP_param.k1*(x1-x2);
    xdot[2] = -MVP_param.k1*x3+MVP_param.k1*MVP_param.SI*x2;
    xdot[3] = -x3*x4+MVP_param.EGP+(x6*MVP_param.km/MVP_param.VG)+MVP_param.k_glu*x9;
    xdot[4] = d-x5*MVP_param.km;
    xdot[5] = (x5-x6)*MVP_param.km;
    xdot[6] = (x4-x7)/MVP_param.tau_GI;
    xdot[7] = u3-x8/MVP_param.tau_glu;
    xdot[8] = (x8-x9)/MVP_param.tau_glu;

    vector<SX> XUD(1,states_arr);
    XUD.push_back(controls_arr);
    XUD.push_back(disturbance_arr);
    f = Function("f", XUD,xdot);

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
    fun_RK4 = Function("F", XUD,vector<SX>{states_arr_RK4});

    SX val_max = log(exp(x1)+exp(x2));
    SX val_min = -log(exp(-x1)+exp(-x2));
    SX val_abs = sqrt(x1*x1+1);
    myMax = Function("myMax",{x1,x2},{val_max});
    myMin = Function("myMin",{x1,x2},{val_min});
    myAbs = Function("myAbs",{x1},{val_abs});

    X = MX::sym("X",n_states*(N+1));
    U = MX::sym("U",n_controls*N);
    P = MX::sym("P",n_states+1+1);
    OPT_variables = vertcat(X,U);


    vector<MX> X_vec(1,X(Slice(0,n_states),Slice()));
    vector<MX> U_vec(1,U(Slice(0,n_controls),Slice()));
    //vector<MX> X_vec(1,X.nz(Slice(0,n_states)));
    //vector<MX> U_vec(1,U.nz(Slice(0,n_controls)));
    for(int i=1;i<N+1;++i){
        X_vec.push_back(X(Slice(i*n_states,(i+1)*n_states),Slice()));
        //X_vec.push_back(X.nz(Slice(i*n_states,(i+1)*n_states)));
        if(i!=N){
            U_vec.push_back(U(Slice(i*n_controls,(i+1)*n_controls),Slice()));
            //U_vec.push_back(U.nz(Slice(i*n_controls,(i+1)*n_controls)));
        }
    }

    obj = 0;
    for(int i=0;i<N;++i){
        MX st = X_vec[i];
        MX zk = st(3), con = U_vec[i];
        MX max_k = vertcat(myMax({zk - z_max, 0}));
        MX min_k = vertcat(myMin({zk - z_min, 0}));
        MX abs_ubo_k = vertcat(myAbs(con(1)));
        MX abs_uglu_k = vertcat(myAbs(con(2)));
        obj = obj + (zk-P(n_states))*Qz*(zk-P(n_states)) + max_k*Qz_max*max_k + min_k*Qz_min*min_k + Qu_bo*abs_ubo_k + Qu_glu*abs_uglu_k;
        obj = obj + (con(0)-uba_ref)*Qu_delta_ba*(con(0)-uba_ref);
        //if(i!=N-1){
        //    MX con_next = U_vec[i+1];
        //    obj = obj + (con_next(0)-con(0))*Qu_delta_ba*(con_next(0)-con(0));
        //}
    }


    vector<MX> g_vec;
    MX st = X_vec[0];
    MX zk = st(3);
    g_vec.push_back(st-P(Slice(0,n_states),Slice()));
    //g_vec.push_back(st-P.nz(Slice(0,n_states)));
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

        g_vec.push_back(st_next-st_next_RK4);
    }
    g = vertcat(g_vec);

    nlp_prob = {{"x", OPT_variables}, {"f", obj}, {"g", g}, {"p",P}};
    //opts["ipopt.tol"] = 1e-8;
    opts["ipopt.max_iter"] = 500;
    opts["ipopt.print_level"] = 3;
    //opts["ipopt.acceptable_tol"] = 1E-8;
    //opts["ipopt.acceptable_obj_change_tol"] = 1E-8;
    solver = nlpsol("nlpsol", "ipopt", nlp_prob, opts);

    vector<double> OPT_varX_min(n_states*(N+1),x_min);
    vector<double> OPT_varX_max(n_states*(N+1),x_max);
    vector<double> OPT_varU_min(n_controls*N,u_min);
    vector<double> OPT_varU_max(n_controls*N,u_max);
    OPT_varX_min.insert(OPT_varX_min.end(),OPT_varU_min.begin(),OPT_varU_min.end());
    OPT_varX_max.insert(OPT_varX_max.end(),OPT_varU_max.begin(),OPT_varU_max.end());

    args["lbx"] = OPT_varX_min;
    args["ubx"] = OPT_varX_max;
    args["lbg"] = vector<double>(n_states*(N+1),0);
    args["ubg"] = vector<double>(n_states*(N+1),0);
}

CasADiSolver::~CasADiSolver() {};
void CasADiSolver::NLPsolverConfig(double *x0, double* xbar,double dmeal, double z_ref) {
    vector<double> OPT_varX0(x0,x0+(N+1)*n_states+N*n_controls);
    vector<double> OPT_p(xbar,xbar+n_states);
    OPT_p.push_back(z_ref);
    OPT_p.push_back(dmeal);
    args["x0"] = OPT_varX0;
    args["p"] = OPT_p;
}
map<string, vector<double>> CasADiSolver::getArgs() {
    return args;
}
Function CasADiSolver::getSolver() {
    return solver;
}