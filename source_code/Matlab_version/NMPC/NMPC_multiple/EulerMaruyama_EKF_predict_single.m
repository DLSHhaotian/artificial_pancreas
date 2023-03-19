function [Xk_RK4,Pk_RK4]=EulerMaruyama_EKF_predict_single(tspan,X,P,G,u,d,p)
% Explicit Euler-Maruyama method with a fixed step size
% Syntax: [Xk_RK4,Pk_RK4]=EulerMaruyama_EKF_predict_single(tspan,X,P,G,u,d,p)
%         tspan: sampling time
%         X: state
%         P: covariance of the state
%         G: matrix for covariance
%         u: input
%         d: disturbance(meal) 
%         p: parameters of MVP model
M = 4; 
DT = tspan/M;
Xk_RK4=X;
Pk_RK4=P;
for i = 1:M
    [X_dot,P_dot] = CDEKFPredict_single_noDiffusion(Xk_RK4, Pk_RK4, u, d, G,p);
    Xk_RK4 = Xk_RK4 + X_dot*DT;
    Xk_RK4(4) = Xk_RK4(4)+p(8)*d;
    Pk_RK4 = Pk_RK4 +P_dot*DT;
end
end
