function [solver,args] = CasADiInit_NMPC_single(N,T,p_MVP,us_MVP,z_max,z_min,z_ref,u_min,u_max,Wz,Wz_max,Wz_min,Wu_bo,Wu_ba)
% Syntax: [solver,args] = CasADiInit_NMPC_single(N,T,p_MVP,us_MVP,z_max,z_min,z_ref,u_min,u_max,Wz,Wz_max,Wz_min,Wu_bo,Wu_ba)
%         N: number of step in control horizon
%         T: sampling time
%         p_MVP: parameters of MVP model
%         us_MVP: steady-state of input of MVP model
%         z_max: upper limit of the controlled variable
%         z_min: lower limit of the controlled variable
%         z_ref: set point of the controlled variable
%         u_min: lower limit of the input
%         u_max: upper limit of the input
%         Wz: weight value of the controlled variable
%         Wz_max: weight value of upper limit of the controlled variable
%         Wz_min: weight value of lower limit of the controlled variable
%         Wu_bo: weight value of the bolus insulin
%         Wu_ba: weight value of the basal insulin

addpath('C:\Users\DLSH\diabate\casadi');
import casadi.*

k1=p_MVP(1);%[min-1] 
CI=p_MVP(2);%[L/min]
SI=p_MVP(3);%[L/mU/min]0.0013
EGP=p_MVP(4);%[mg/dL/min]
km=p_MVP(5);%[min-1] meal absorption tau_D 1/tau_G
VG=p_MVP(6);%[dL] VG[L]
tau_GI=p_MVP(7);%[min]
rou4=p_MVP(8);%[mg/dL/min]

% Declare model variables
x1=SX.sym('x1');
x2=SX.sym('x2');
x3=SX.sym('x3');
x4=SX.sym('x4');
x5=SX.sym('x5');
x6=SX.sym('x6');
x7=SX.sym('x7');
x8=SX.sym('x8');
x9=SX.sym('x9');

states=[x1;x2;x3;x4;x5;x6;x7];
n_states=length(states);

u1=SX.sym('u1');%basal
u2=SX.sym('u2');%bolus
controls=[u1;u2]; 
n_controls=length(controls);

d=SX.sym('d');%meal intake
disturbance=d;

rhs=[k1*((u1+u2)/CI-x1);
    k1*(x1-x2);
    -k1*x3+k1*SI*x2;
    -x3*x4+EGP+(x6*km/VG);
    d-x5*km;
    (x5-x6)*km;
    (x4-x7)/tau_GI]; 

x_max=log(exp(x1)+exp(x2));
x_min=-log(exp(-x1)+exp(-x2));
x_abs=sqrt(x1*x1+1);
myMin=Function('myMin',{x1,x2},{x_min});
myMax=Function('myMax',{x1,x2},{x_max});
myAbs=Function('myAbs',{x1},{x_abs});

f=Function('f',{states,controls,disturbance},{rhs}); % nonlinear mapping function f(x,u)
U=SX.sym('U',n_controls,N); % Decision variables (controls)
P=SX.sym('P',n_states + 1 + 1);%% parameters [x0;z_ref;d](which include the initial state and the reference output and disturbance)
X=SX.sym('X',n_states,(N+1));

% Fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = T/M;
states_RK4=states;
for j=1:M
    RK1 = f(states_RK4, controls,disturbance);
    RK2 = f(states_RK4 + DT/2 * RK1, controls,disturbance);
    RK3 = f(states_RK4 + DT/2 * RK2, controls,disturbance);
    RK4 = f(states_RK4 + DT * RK3, controls,disturbance);
    states_RK4=states_RK4+DT/6*(RK1 +2*RK2 +2*RK3 +RK4);
end
fun_RK4 = Function('F', {states, controls, disturbance}, {states_RK4}, {'x0','u','d'}, {'xf'});

obj=0; % Objective function
g=[]; % constraints vector
Qz=Wz; % weighing matrices (output z)
Qz_max=Wz_max;% weighing matrices (output z penalty functions for max)
Qz_min=Wz_min;% weighing matrices (output z penalty functions for min)
R_bo=Wu_bo;% weighing matrices (controls bolus)
R_ba=Wu_ba;% weighing matrices (controls basal)
uba_ref = us_MVP(1);

st=X(:,1); % initial state
zk=st(4); % output z
g=[g;st-P(1:n_states)]; % initial condition constraints，

for k=1:N
    st=X(:,k); zk=st(4);
    con=U(:,k); 
    obj = obj+(zk-P(n_states+1))*Qz*(zk-P(n_states+1)) + myMax(zk-z_max,0)*Qz_max*myMax(zk-z_max,0) + myMin(zk-z_min,0)*Qz_min*myMin(zk-z_min,0)+R_bo*myAbs(con(2));
    if k~=N
        %delta uba 2-norm
%          con_next=U(:,k+1);
%          obj=obj+(con_next(1)-con(1))*R_ba*(con_next(1)-con(1));
    end
    %uba-uba_bar 2-norm
    obj = obj+(con(1)-uba_ref)*R_ba*(con(1)-uba_ref);
    st_next=X(:,k+1);
    if k==1    
        obj = obj+(con(1)-uba_ref)*R_ba*(con(1)-uba_ref);
        st_next_RK4=fun_RK4(st,con,P(n_states+2));
    else
        st_next_RK4=fun_RK4(st,con,0);
    end
    g=[g;st_next-st_next_RK4]; % compute constraints
end

OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
opts = struct;
opts.ipopt.max_iter = 500;
opts.ipopt.print_level =3;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-8;
solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:n_states*(N+1),1) = 0; % -1e-20 % Equality constraints  
args.ubg(1:n_states*(N+1),1) = 0; % 1e-20 % Equality constraints
args.lbx(1:n_states*(N+1),1) = -inf; %state x lower bound
args.ubx(1:n_states*(N+1),1) = inf; %state x upper bound
args.lbx(n_states*(N+1)+1:n_states*(N+1)+n_controls*N,1) = u_min; %v lower bound
args.ubx(n_states*(N+1)+1:n_states*(N+1)+n_controls*N,1) = u_max; %v upper bound

end


