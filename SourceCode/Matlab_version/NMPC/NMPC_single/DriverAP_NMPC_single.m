% Matlab version
% Driver file to run the closed-loop simulation test for artificial pancreas(AP) 
% Nonlinear model predictive control(LMPC) for single-hormone 

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
addSensorNoise = false;
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
u0 = [0;0]; 
ds = 0;
ys=108/18;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_Hovorka = fsolve(@HovorkaModelSteadyState_single,u0,options,ys,ds,p_Hovorka)
xs_Hovorka = HovorkaModely2x_single(us_Hovorka,ys,ds,p_Hovorka)
% --------------------------------------------------------------
% Steady State of extended MVP model(control model)
% --------------------------------------------------------------
u0 = [0;0]; 
ds = 0; 
ys=108;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_MVP = fsolve(@MVPModelSteadyState_single,u0,options,ys,ds,p_MVP)
xs_MVP = MVPModely2x_single(us_MVP,ys,ds,p_MVP)

%%
% --------------------------------------------------------------
% CasADi Initialization
% --------------------------------------------------------------

%glucose concentration
z_min=70; %[mg/dL] hypoglycemia (BG<70 mg/dL) (hypoglycemia is highly undesirable)
z_max=180; %[mg/dL] hyperglycemia (BG>180 mg/dL)
z_ref=108; %[mg/dL] glucose setpoint 

%basal and bolus insulin infusion
u_min=0; %[mU/min]
u_max=inf; %[mU/min]

%state
x_min = -inf;
x_max = inf;

%weight of object function
Wz=0.01;
Wz_max=1;
Wz_min=300000;
Wu_bo=100;
Wu_ba=150;

% NLP solver Initialization
[solver,args] = CasADiInit_NMPC_single(N,T,p_MVP,us_MVP,z_max,z_min,z_ref,u_min,u_max,Wz,Wz_max,Wz_min,Wu_bo,Wu_ba);

%%
% --------------------------------------------------------------
% Standard Wiener process(process noise)
% --------------------------------------------------------------
Ns = 11; % Number of realizations
seed = 100; % Seed for reproducibility
[W,Tw,dW] = ScalarStdWienerProcess(t_final,N_all,Ns,seed);

%%
% --------------------------------------------------------------
% Closed-loop simulation(SDE)
% --------------------------------------------------------------
[x_arr,z_arr,y_arr,u_arr,d_arr,Z_bar,t_cost] = ClosedLoopSimulation_NMPC_single(xs_Hovorka,xs_MVP,us_MVP,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,solver,args);

%%
% --------------------------------------------------------------
% Plot the result 
% --------------------------------------------------------------
PlotAP_single(x_arr,z_arr,y_arr,u_arr,d_arr,Z_bar,t_cost);
