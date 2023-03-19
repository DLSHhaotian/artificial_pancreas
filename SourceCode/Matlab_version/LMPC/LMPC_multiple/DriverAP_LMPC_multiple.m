% Matlab version
% Driver file to run the closed-loop simulation test for artificial pancreas(AP) 
% Linear model predictive control(LMPC) for multiple-hormone 

%%
clc
clear
close all

%%
% --------------------------------------------------------------
% Initialization of closed-loop simulation
% --------------------------------------------------------------
t0=0; % Starting time
T=5; % [min] Sampling time
t_final=24*60; % End time
t = [t0:T:t_final]'; % Sample instants
N_all = length(t);
N = 72; % Steps of control horizon

% --------------------------------------------------------------
% ODE: addProcessNoise = false
% SDE: addProcessNoise = true
% --------------------------------------------------------------
addSensorNoise = true;
addProcessNoise = false;

%%
% --------------------------------------------------------------
% Parameter Estimation
% --------------------------------------------------------------
numPatient = 1;
[t_arr,d_arr,u_arr,y_arr,p_Hovorka,p_MVP,xs_MVP,dW] = ParamEstimate_dataset(numPatient);
p_MVP = ParamEstimateMLE(t_arr,u_arr,d_arr,y_arr,dW,T,p_MVP,xs_MVP);

%%
% --------------------------------------------------------------
% Steady State of extended Hovorka model(simulation model)
% --------------------------------------------------------------
u0 = [0;0;0]; 
ds = 0;
HRs = 70;
ys = 108/18;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_Hovorka = fsolve(@HovorkaModelSteadyState_multiple,u0,options,ys,ds,p_Hovorka,HRs)
xs_Hovorka = HovorkaModely2x_multiple(us_Hovorka,ys,ds,p_Hovorka)
% --------------------------------------------------------------
% Steady State of extended MVP model(control model)
% --------------------------------------------------------------
u0 = [0;0;0]; 
ds = 0; 
ys = 108;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_MVP = fsolve(@MVPModelSteadyState_multiple,u0,options,ys,ds,p_MVP)
xs_MVP = MVPModely2x_multiple(us_MVP,ys,ds,p_MVP)

%%
% --------------------------------------------------------------
% Linearization and discretization of MVP model of single-hormone
% --------------------------------------------------------------
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 u1 u2 u3 d1 
X = [x1;x2;x3;x4;x5;x6;x7;x8;x9]; 
F = [u1;u2;u3;d1]; 
%p_MVP=[k1_MVP; CI_MVP; SI_MVP; EGP_MVP; km_MVP; VG_MVP; tau_GI_MVP; rou4_MVP];
tau_glu = 19;
k_glu = 0.075;

fx(1,1) = p_MVP(1)*((u1+u2)/p_MVP(2)-x1);
fx(2,1) = p_MVP(1)*(x1-x2);
fx(3,1) = -p_MVP(1)*x3+p_MVP(1)*p_MVP(3)*x2;
fx(4,1) = -x3*x4+p_MVP(4)+(x6*p_MVP(5)/p_MVP(6))+(k_glu*x9);
fx(5,1) = d1-x5*p_MVP(5);
fx(6,1) = (x5-x6)*p_MVP(5);
fx(7,1) = (x4-x7)/p_MVP(7);
fx(8,1) = u3-x8/tau_glu;
fx(9,1) = (x8-x9)/tau_glu;
gx = x4;

Ass=eval(vpa(subs(jacobian(fx,X),[X;F],[xs_MVP;us_MVP;ds]),4));
Bss=eval(vpa(subs(jacobian(fx,F(1:3)),[X;F],[xs_MVP;us_MVP;ds]),4));
Bdss=eval(vpa(subs(jacobian(fx,F(4)),[X;F],[xs_MVP;us_MVP;ds]),4));
Css=eval(vpa(subs(jacobian(gx,X),[X;F],[xs_MVP;us_MVP;ds]),4));

%Compute discrete-time state space models
[Ad,Bd]=c2dzoh(Ass,Bss,T)
[Ad,Bd_d]=c2dzoh(Ass,Bdss,T)

%%
% --------------------------------------------------------------
% CasADi Initialization
% --------------------------------------------------------------

%glucose concentration
z_min=80; %[mg/dL] hypoglycemia (BG<70 mg/dL) (hypoglycemia is highly undesirable)
z_max=180; %[mg/dL] hyperglycemia (BG>180 mg/dL)
z_ref=108; %[mg/dL] glucose setpoint 

%basal and bolus insulin infusion
u_min=0; %[mU/min]
u_max=inf; %[mU/min]

%state
x_min = -inf;
x_max = inf;

%weight of object function
Wz = 0.01;
Weta_min = 2; %z_min
Weta_max = 10; %z_max
Wu_ba = 100;
Wu_delta_ba = 0;
Wu_bo = 100;
Wu_glu = 12000;

solver = CasADiInit_LMPC_multiple(N,T,Ad,Bd,Bd_d,Css,xs_MVP,us_MVP,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo,Wu_glu);

%%
% --------------------------------------------------------------
% Standard Wiener process(process noise)
% --------------------------------------------------------------
Ns = 16; % Number of realizations
seed = 100; % Seed for reproducibility
[W,Tw,dW] = ScalarStdWienerProcess(t_final,N_all,Ns,seed);

%%
% --------------------------------------------------------------
% Closed-loop simulation
% --------------------------------------------------------------
[x_arr,z_arr,y_arr,u_arr,d_arr,HR_arr,Z_bar,t_cost] = ClosedLoopSimulation_LMPC_multiple(Ad,Bd,Bd_d,Css,solver,xs_Hovorka,xs_MVP,us_MVP,ys,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo,Wu_glu);

%%
% --------------------------------------------------------------
% Plot the result 
% --------------------------------------------------------------
PlotAP_multiple(x_arr,z_arr,y_arr,u_arr,d_arr,HR_arr,Z_bar,t_cost);
