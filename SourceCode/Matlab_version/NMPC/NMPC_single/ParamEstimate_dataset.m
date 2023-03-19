function [t_arr,d_arr,u_arr,y_arr,p_Hovorka,p_MVP,xs_MVP,dw] = ParamEstimate_dataset(numPatient)
% Syntax: [t_arr,d_arr,u_arr,y_arr,p_Hovorka,p_MVP,xs_MVP,dw] = ParamEstimate_dataset(numPatient)
%         numPatient: number of patient

[seed,num] = HovorkaPatientSeed(numPatient);

% Generate virtual patient parameters
p_Hovorka = HovorkaParam_random(seed);

% --------------------------------------------------------------
% Steady State of extended Hovorka model(simulation model)
% --------------------------------------------------------------
u0 = [0;0]; 
ds = 0;
ys=108/18;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_Hovorka = fsolve(@HovorkaModelSteadyState_single,u0,options,ys,ds,p_Hovorka);
xs_Hovorka = HovorkaModely2x_single(us_Hovorka,ys,ds,p_Hovorka);


% --------------------------------------------------------------
% Initial guess of MVP model parameters
% --------------------------------------------------------------
k1_MVP=0.0182;%[min-1] 
CI_MVP=1.5;%[L/min]
SI_MVP=0.00092;%[L/mU/min]
EGP_MVP=1.25;%[mg/dL/min]
km_MVP=0.027;%[min-1] 
VG_MVP=200;%[dL] 
tau_GI_MVP=6.7;%[min]
rou4_MVP=3.05;%[mg/dL/min]
p_MVP=[k1_MVP; CI_MVP; SI_MVP; EGP_MVP; km_MVP; VG_MVP; tau_GI_MVP; rou4_MVP];

% --------------------------------------------------------------
% Steady State of extended MVP model(control model)
% --------------------------------------------------------------
u0 = [0;0]; 
ds = 0; 
ys=108;
options=optimoptions("fsolve",'Algorithm','levenberg-marquardt');
us_MVP = fsolve(@MVPModelSteadyState_single,u0,options,ys,ds,p_MVP);
xs_MVP = MVPModely2x_single(us_MVP,ys,ds,p_MVP);

% --------------------------------------------------------------
% Initialization of closed-loop simulation
% --------------------------------------------------------------
t0=0; % Starting time
T=5; % [min] Sampling time
N =72;
t_final=30*60; % End time
t_arr = [t0:T:t_final]'; % Sample instants
N_all = length(t_arr);

% --------------------------------------------------------------
% Standard Wiener process(process noise)
% --------------------------------------------------------------
Ns = 11; % Number of realizations
seed = 100; % Seed for reproducibility
[W,Tw,dW] = ScalarStdWienerProcess(t_final,N_all,Ns,seed);
dw = dW(num,:);

% Number of states, inputs, etc. (Hovorka model) 
nx = 11; nu = 2; ny = 1; nz = 1; nd = 1;
% Number of states, inputs, etc. (MVP model)
n_states = 7; n_controls = 2; n_disturbance = 1;
% Set up inputs
u_arr = zeros(nu,N_all-1);
for i = 1:10
    u_arr(1,(i-1)*36+1:i*36) = us_Hovorka(1)+normrnd(0,sqrt(0.3*us_Hovorka(1)));
end
ubo_test = 1000;
u_arr(2,8*60/T)=ubo_test; u_arr(2,11.5*60/T)=ubo_test; u_arr(2,13.25*60/T)=ubo_test; u_arr(2,18*60/T)=ubo_test; u_arr(2,22*60/T)=ubo_test;
% Set up disturbance(meals)
d_arr = zeros(nd,N_all-1);
d_arr(:,8*60/T)=72/T; d_arr(:,11.5*60/T)=36/T; d_arr(:,13.25*60/T)=131/T; d_arr(:,18*60/T)=51/T; d_arr(:,22*60/T)=70/T;

% Set up cost(calculation time)
t_cost=zeros(1,N_all-1);

% Set up output and states
y_arr = zeros(ny,N_all);
z_arr = zeros(nz,N_all);
x0 = xs_Hovorka;
x_arr = zeros(nx,N_all);
x_arr(:,1) = x0;


% --------------------------------------------------------------
% Measurement noise
% --------------------------------------------------------------
noise_cc_prepre=0;
noise_cc_pre=0;
noise_v_prepre=0;
noise_v_pre=0;

R_wcc = 11.3;
Lr_wcc = chol(R_wcc,'lower');
wcc_arr = Lr_wcc*randn(1,N_all);
R_w = 14.45;
Lr_w = chol(R_w,'lower');
w_arr = Lr_w*randn(1,N_all);

% --------------------------------------------------------------
% Closed-loop simulation(SDE)
% --------------------------------------------------------------
for k = 1:N_all-1
% Update of measurement noise    
[noise_cc,noise_v] = sensorNoise(noise_cc_prepre,noise_cc_pre,noise_v_prepre,noise_v_pre,wcc_arr(k),w_arr(k));
noise_cc_prepre=noise_cc_pre;
noise_cc_pre=noise_cc;
noise_v_prepre=noise_v_pre;
noise_v_pre=noise_v;

% Update of output
y_arr(:,k) = HovarkaSensor(x_arr(:,k),noise_cc+noise_v);
z_arr(:,k) = HovarkaOutput(x_arr(:,k)); 

% Simulation of Hovorka model with single hormone
uk=u_arr(:,k);
xk = EulerMaruyamaSDE_Hovorka_single(T,1.25,x_arr(:,k),[uk(1);uk(2)],d_arr(:,k),dW(:,k)./T,xs_Hovorka*0.01,p_Hovorka);
x_arr(:,k+1) = xk;
end
k=N_all;
y_arr(:,k) = HovarkaSensor(x_arr(:,k),noise_cc+noise_v);
z_arr(:,k) = HovarkaOutput(x_arr(:,k));

end