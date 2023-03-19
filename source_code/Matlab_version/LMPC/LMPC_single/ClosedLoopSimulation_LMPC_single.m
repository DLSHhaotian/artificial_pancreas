function [x_arr,z_arr,y_arr,u_arr,d_arr,Z_bar,t_cost] = ClosedLoopSimulation_LMPC_single(Ad,Bd,Bd_d,Css,solver,xs_Hovorka,xs_MVP,us_MVP,ys,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo)
% Syntax: [x_arr,z_arr,y_arr,u_arr,d_arr,Z_bar,t_cost] = ClosedLoopSimulation_LMPC_single(Ad,Bd,Bd_d,Css,solver,xs_Hovorka,xs_MVP,us_MVP,ys,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo)
%         Ad: state matrix of MVP model
%         Bd: input matrix of MVP model
%         Bd_d: disturbance matrix of MVP model
%         Css: output matrix of MVP model
%         solver: solver of OCP by CasADi
%         xs_Hovorka: steady-state of state of Hovorka model
%         xs_MVP: steady-state of state of MVP model
%         us_MVP: steady-state of input of MVP model
%         N: number of step in control horizon
%         N_all: number of step in simulation
%         T: sampling time
%         t: simulation time vector
%         p_Hovorka: parameters of Hovorka model
%         p_MVP: parameters of MVP model
%         addSensorNoise: add the sensor noise or not
%         dW: process noise
%         addProcessNoise: add the process noise or not
%         z_min: steady-state of input of MVP model
%         z_max: number of step in control horizon
%         z_ref: number of step in simulation
%         u_min: sampling time
%         u_max: simulation time vector
%         Wz: weight value of the controlled variable
%         Weta_min: weight value of lower limit of the controlled variable(slack variable)
%         Weta_max: weight value of upper limit of the controlled variable(slack variable)
%         Wu_ba: weight value of the basal insulin
%         Wu_delta_ba: weight value of the change rate of basal insulin
%         Wu_bo: weight value of the bolus insulin


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

% Number of states, inputs, etc. (Hovorka model) 
nx = 11; nu = 2; ny = 1; nz = 1; nd = 1;
% Number of states, inputs, etc. (MVP model)
n_states = 7; n_controls = 2; n_disturbance = 1;
% Set up inputs
u_arr = zeros(nu,N_all-1);

% Set up disturbance(meals)
d_arr = zeros(nd,N_all-1);
d_arr(:,6*60/T)=50/T; d_arr(:,12*60/T)=60/T; d_arr(:,18*60/T)=70/T;

% Set up cost(calculation time)
t_cost=zeros(1,N_all-1);

% Set up output and states
y_arr = zeros(ny,N_all);
z_arr = zeros(nz,N_all);
x0 = xs_Hovorka;
x_arr = zeros(nx,N_all);
x_arr(:,1) = x0;

% --------------------------------------------------------------
% Initialization of CDEKF
% --------------------------------------------------------------
G=zeros(n_states,1);
G(4)=p_MVP(8);
P_kf = G*G';
R_kf = 0.02^2;
W_kf=diag([1 1 1 0.001 1 1 1])*T;
G_w=diag([1 1 1 1 1 1 1]);


Z_bar = zeros(1,N_all);
y_bar = zeros(ny,N_all);
X_bar = zeros(n_states,N_all);
X_bar(:,1) = xs_MVP;
x_bar = X_bar(:,1);


% --------------------------------------------------------------
% Closed-loop simulation(SDE)
% --------------------------------------------------------------
% Initialization of the input decision variables
u0 = zeros(N,n_controls); 
% Initialization of the slack decision variables
eta0 = zeros(2*N,1);
uk = us_MVP;
h = waitbar(0,'Please wait...');
for k = 1:N_all-1
% Update of measurement noise    
[noise_cc,noise_v] = sensorNoise(noise_cc_prepre,noise_cc_pre,noise_v_prepre,noise_v_pre,wcc_arr(k),w_arr(k));
noise_cc_prepre=noise_cc_pre;
noise_cc_pre=noise_cc;
noise_v_prepre=noise_v_pre;
noise_v_pre=noise_v;

% Update of output
if(addSensorNoise)
    y_arr(:,k) = HovarkaSensor(x_arr(:,k),noise_cc+noise_v);
else
    y_arr(:,k) = HovarkaSensor(x_arr(:,k),0);
end
z_arr(:,k) = HovarkaOutput(x_arr(:,k)); 

tic;
% KF filting and prediction
x_bar=X_bar(:,k);
[x_bar,z_bar,P_kf] = KF_single(T,uk,x_bar,P_kf,R_kf,W_kf,G_w,d_arr(:,k)*1000,y_arr(:,k),Ad,Bd,Bd_d,Css,us_MVP,xs_MVP,ys);
X_bar(:,k+1)=x_bar+xs_MVP;
Z_bar(:,k)=z_bar+ys;

% Solve optimal control problem
% single-shooting
OPT_var_x0 = [reshape(u0',n_controls*N,1);eta0];
sol = CasADiCompute_LMPC_single(N,T,Ad,Bd,Bd_d,Css,xs_MVP,us_MVP,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo,x_bar-xs_MVP,d_arr(k)*1000,uk-us_MVP,OPT_var_x0,solver);

% Only get controls from the solution
u_hor = reshape(sol(1:N*n_controls)',n_controls,N);
% Compute the cost(calculation time)
elapsedTime = toc;
t_cost(:,k)=elapsedTime;
u_arr(1,k)=u_hor(1,1)+us_MVP(1);
u_arr(2,k)=u_hor(2,1)+us_MVP(2);
% Construction of input decision variables for next time point
u0 = [u_hor(:,2:size(u_hor,2)),u_hor(:,size(u_hor,2))]';
% Construction of slack decision variables for next time point
eta0 = sol(N*n_controls+1:end);


% Simulation of Hovorka model with single hormone
uk=u_arr(:,k);
if(addProcessNoise)
    xk = EulerMaruyamaSDE_Hovorka_single(T,1.25,x_arr(:,k),[uk(1);uk(2)],d_arr(:,k),dW(:,k)./T,xs_Hovorka*0.01,p_Hovorka);
else
    [Tk,Xk] = ode45(@HovorkaModel_single,[t(k) t(k+1)],x_arr(:,k),[],[uk(1);uk(2)],d_arr(:,k),p_Hovorka);
    xk=Xk(end,:)';
end
x_arr(:,k+1) = xk;

per = k / (N_all-1);
waitbar(per, h ,sprintf('%2.0f%%',per*100))
end
k=N_all;
y_arr(:,k) = HovarkaSensor(x_arr(:,k),noise_cc+noise_v);
z_arr(:,k) = HovarkaOutput(x_arr(:,k));
Z_bar(:,k)=z_bar;
close(h)
end