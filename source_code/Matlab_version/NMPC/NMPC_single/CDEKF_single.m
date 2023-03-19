function [x_bar,z_bar,P_kf] = CDEKF_single(T,uk,x_bar,P_kf,R_kf,G,dk,yk,p_MVP)
% Syntax: [x_bar,z_bar,P_kf] = CDEKF_single(T,uk,x_bar,P_kf,R_kf,G,dk,yk,p_MVP)
%         T: sampling time
%         uk: current calculated input 
%         x_bar: observations of the state
%         P_kf: covariance of the state
%         R_kf: covariance of the error
%         G: matrix for covariance
%         dk: current disturbance(meal intake)
%         yk: output measured by CGM
%         p_MVP: parameters of MVP model

% --------------------------------------------------------------
% Prediction
% --------------------------------------------------------------
[x_bar,P_kf]=RK4_EKF_predict_single(T,x_bar,P_kf,G,[uk;dk],0,p_MVP);%MVP model meal[mg/min]
% --------------------------------------------------------------
% Update
% --------------------------------------------------------------
[C_kf,z_bar] = MVPOutputJacobian_single(x_bar);
e_kf=yk-z_bar;
R_kf=C_kf*P_kf*C_kf'+R_kf;
K_kf=P_kf*C_kf'/R_kf;
x_bar=x_bar+K_kf*e_kf;
P_kf=P_kf-K_kf*R_kf*K_kf';
end
