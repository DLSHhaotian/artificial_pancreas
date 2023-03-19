function [x_arr,z_arr,y_arr,u_arr,d_arr,HR_arr,Z_bar,t_cost] = ClosedLoopSimulation_NMPC_multiple(xs_Hovorka,xs_MVP,us_MVP,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,solver,args)
% Syntax: [x_arr,z_arr,y_arr,u_arr,d_arr,HR_arr,Z_bar,t_cost] = ClosedLoopSimulation_NMPC_multiple(xs_Hovorka,xs_MVP,us_MVP,N,N_all,T,t,p_Hovorka,p_MVP,addSensorNoise,dW,addProcessNoise,solver,args)
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
%         solver: solver of OCP by CasADi
%         args: info of OCP 
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
nx = 16; nu = 3; ny = 1; nz = 1; nd = 1;
% Number of states, inputs, etc. (MVP model)
n_states = 9; n_controls = 3; n_disturbance = 1;
% Set up inputs
u_arr = zeros(nu,N_all-1);

% Set up disturbance(meals)
d_arr = zeros(nd,N_all-1);
d_arr(:,6*60/T)=50/T; d_arr(:,12*60/T)=60/T; d_arr(:,18*60/T)=70/T;

% Set up exercise
HR_base = 70;
HR_arr = HR_base*ones(1,N_all-1);
HR_arr(:,15*60/T:15.75*60/T)=150;

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
% Initialization of the states decision variables
X0 = repmat(x_bar,N+1,1); 
uk = zeros(3,1);
z_ref = 108;
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
% CDEKF filting and prediction
x_bar=X_bar(:,k);
[x_bar,z_bar,P_kf] = CDEKF_multiple(T,uk,x_bar,P_kf,R_kf,G,d_arr(:,k)*1000,y_arr(:,k),p_MVP);
X_bar(:,k+1)=x_bar;
Z_bar(:,k)=z_bar;

% Solve optimal control problem
% Set the values of the parameters vector
args.p = [x_bar;z_ref;d_arr(k)*1000];
% Set the initial value of the decision variables
args.x0 = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1)];
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
% Only get controls from the solution
u_hor = reshape(full(sol.x(n_states*(N+1)+1:end))',n_controls,N);
% Compute the cost(calculation time)
elapsedTime = toc;
t_cost(:,k)=elapsedTime;

u_arr(:,k)=u_hor(:,1);
% Construction of input decision variables for next time point
u0 = [u_hor(:,2:size(u_hor,2)),u_hor(:,size(u_hor,2))]';
% Construction of state decision variables for next time point
X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; 
X0 = [X0(2:end,:);X0(end,:)];


% Simulation of Hovorka model with multiple hormone
uk=u_arr(:,k);
if(addProcessNoise)
    xk = EulerMaruyamaSDE_Hovorka_multiple(T,1.25,x_arr(:,k),[uk(1);uk(2);uk(3)],d_arr(:,k),dW(:,k)./T,xs_Hovorka*0.01,p_Hovorka,HR_arr(:,k));
else
    [Tk,Xk] = ode45(@HovorkaModel_multiple,[t(k) t(k+1)],x_arr(:,k),[],[uk(1);uk(2);uk(3)],d_arr(:,k),p_Hovorka,HR_arr(:,k));
    xk=Xk(end,:)';
end
x_arr(:,k+1) = xk;

per = k / (N_all-1);
waitbar(per, h ,sprintf('%2.0f%%',per*100))
end
k=N_all;
y_arr(:,k) = HovarkaSensor(x_arr(:,k),noise_cc+noise_v);
z_arr(:,k) = HovarkaOutput(x_arr(:,k));

close(h)
end