function [x_bar,z_bar,P_kf] = KF_single(T,uk,x_bar,P_kf,R_kf,W_kf,G_w,dk,yk,Ad,Bd,Bd_d,Css,us_MVP,xs_MVP,ys)
% Syntax: [x_bar,z_bar,P_kf] = KF_single(T,uk,x_bar,P_kf,R_kf,W_kf,G_w,dk,yk,Ad,Bd,Bd_d,Css,us_MVP,xs_MVP,ys)
%         T: sampling time
%         uk: current calculated input 
%         x_bar: observations of the state
%         P_kf: covariance of the state
%         R_kf: covariance of the error
%         W_kf: matrix for covariance
%         G_w: gain matrix for covariance
%         dk: current disturbance(meal intake)
%         yk: output measured by CGM
%         Ad: state matrix of MVP model
%         Bd: input matrix of MVP model
%         Bd_d: disturbance matrix of MVP model
%         Css: output matrix of MVP model
%         us_MVP: steady-state of input of MVP model
%         xs_MVP: steady-state of state of MVP model
%         ys: steady-state of output of MVP model

uk_tmp = zeros(2,1);
uk_tmp(1) = uk(1)-us_MVP(1);
uk_tmp(2) = uk(2)-us_MVP(2);
% --------------------------------------------------------------
% Prediction
% --------------------------------------------------------------
x_bar=x_bar-xs_MVP;
x_bar=Ad*x_bar+Bd*uk_tmp+Bd_d*dk;
P_kf=Ad*P_kf*Ad'+G_w*W_kf*G_w';
% --------------------------------------------------------------
% Update
% --------------------------------------------------------------
z_bar=Css*x_bar;
e_kf=yk-z_bar-ys;
R_kf=Css*P_kf*Css'+R_kf;
K_kf=P_kf*Css'/R_kf;
x_bar=x_bar+K_kf*e_kf;
P_kf=P_kf-K_kf*R_kf*K_kf';
end
