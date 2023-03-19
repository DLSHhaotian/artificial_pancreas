function u_opt = CasADiCompute_LMPC_single(N,T,Ad,Bd,Bd_d,Css,xs_MVP,us_MVP,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo,x_bar,dmeal,u_1,OPT_var_x0,solver)
% Syntax: u_opt = CasADiCompute_LMPC_single(N,T,Ad,Bd,Bd_d,Css,xs_MVP,us_MVP,z_min,z_max,z_ref,u_min,u_max,Wz,Weta_min,Weta_max,Wu_ba,Wu_delta_ba,Wu_bo,x_bar,dmeal,u_1,OPT_var_x0,solver)
%         N: number of step in control horizon
%         T: sampling time
%         Ad: state matrix of MVP model
%         Bd: input matrix of MVP model
%         Bd_d: disturbance matrix of MVP model
%         Css: output matrix of MVP model
%         xs_MVP: steady-state of state of MVP model
%         us_MVP: steady-state of input of MVP model
%         z_min: lower limit of the controlled variable
%         z_max: upper limit of the controlled variable
%         z_ref: set point of the controlled variable
%         u_min: lower limit of the input
%         u_max: upper limit of the input
%         Wz: weight value of the controlled variable
%         Weta_min: weight value of lower limit of the controlled variable(slack variable)
%         Weta_max: weight value of upper limit of the controlled variable(slack variable)
%         Wu_ba: weight value of the basal insulin
%         Wu_delta_ba: weight value of the change rate of basal insulin
%         Wu_bo: weight value of the bolus insulin
%         x_bar: observation of state
%         dmeal: meal intake
%         u_1: calculated input last time
%         OPT_var_x0: initial guess of decision variable
%         solver: solver of OCP by CasADi

addpath('C:\Users\DLSH\diabate\casadi');
import casadi.*


num_x=size(Ad,1);
num_u=size(Bd,2);
num_d=size(Bd_d,2);
num_y=1;
num_z=1;
num_eta=num_y;
num_u_delta_cons=N;
num_z_cons=N*num_y;

Phi=zeros(N*num_y,num_x);
Phi_d=zeros(N*num_y,num_d);
Gam=zeros(N*num_y,N*num_u+2*N*num_eta);
%Gam_d=zeros(N*num_y,N*num_d);
%Weight matrix for tracking the reference Qz
Qz=diag(Wz*ones(N*num_y,1));
%Weight matrix for input rate S
Sz=diag([Wu_delta_ba;0]);
%Quadratic objective term Hs
Hs=zeros(N*num_u+2*N*num_eta,N*num_u+2*N*num_eta);
%Linear objective term Mu_delta
Mu_delta=zeros(N*num_u+2*N*num_eta,num_u);
Mu_delta(1:num_u,:)=-Sz;
%Quadratic objective term H_eta
H_eta1=diag([zeros(N*num_u,1);Weta_min*ones(N*num_eta,1);zeros(N*num_eta,1)]);
H_eta2=diag([zeros(N*num_u,1);zeros(N*num_eta,1);Weta_max*ones(N*num_eta,1)]);

H_ba = diag(zeros(N*num_u+2*N*num_eta,1));
H_bo = diag(zeros(N*num_u+2*N*num_eta,1));
gbo = zeros(N*num_u+2*N*num_eta,1);
for i=1:N
    %fill the Phi
    Phi((i-1)*num_y+1:i*num_y,:)=Css*Ad^i;
    Phi_d((i-1)*num_y+1:i*num_y,:)=Css*Ad^i*Bd_d;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam
    Gam_i=zeros(num_y,N*num_u+2*N*num_eta);%[num_x * num_u*N]
    for j=1:i
        Gam_i(:,(j-1)*num_u+1:j*num_u)=Css*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Gam((i-1)*num_y+1:i*num_y,:)=Gam_i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:i*num_u)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_ba((i-1)*num_u+1,(i-1)*num_u+1) = Wu_ba;
    gbo((i)*num_u,1) = Wu_bo;
    H_bo(i*num_u,i*num_u) = Wu_bo;
end
Mx0=Gam'*Qz*Phi;
Md=Gam'*Qz*Phi_d;
%Linear objective term Mr
Mr=-Gam'*Qz;
%Quadratic objective term Hr
Hr=Gam'*Qz*Gam;
Hu=Hr+Hs+H_eta1+H_eta2+H_ba;
g=Mx0*x_bar+Mr*(z_ref*ones(N,1))+Mu_delta*u_1+Md*dmeal+gbo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Az_cons_min=Gam;
Az_cons_min(:,N*num_u+1:N*num_u+N*num_eta)=diag(1*ones(N*num_eta,1));
Az_cons_min(:,N*num_u+N*num_eta+1:end)=diag(zeros(N*num_eta,1));
Az_cons_max=Gam;
Az_cons_max(:,N*num_u+1:N*num_u+N*num_eta)=diag(zeros(N*num_eta,1));
Az_cons_max(:,N*num_u+N*num_eta+1:end)=diag(-1*ones(N*num_eta,1));
Z_min=repmat(z_min,N,1)-Phi*x_bar-Phi_d*dmeal;
Z_max_fake=5000*ones(num_z_cons,1);
Z_max=repmat(z_max,N,1)-Phi*x_bar-Phi_d*dmeal;
Z_min_fake=-5000*ones(num_z_cons,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPT_var_min=[u_min*ones(N*num_u,1);zeros(2*N*num_eta,1)];
OPT_var_max=[u_max*ones(N*num_u,1);inf*ones(2*N*num_eta,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_ieq=[Az_cons_min;Az_cons_max];
b_ieq_min=[Z_min;Z_min_fake];
b_ieq_max=[Z_max_fake;Z_max];

Hu_DM = DM(Hu);
g_DM = DM(g);
OPT_var_min_DM = DM(OPT_var_min);
OPT_var_max_DM = DM(OPT_var_max);
A_ieq_DM = DM(A_ieq);
b_ieq_min_DM = DM(b_ieq_min);
b_ieq_max_DM = DM(b_ieq_max);

r = solver('h',Hu_DM,'g',g_DM,'a',A_ieq_DM,'lbx',OPT_var_min_DM,'ubx',OPT_var_max_DM,'lba', b_ieq_min_DM,'uba',b_ieq_max_DM);
u_opt = full(r.x);
end