function [Xk_RK4,Pk_RK4]=RK4_EKF_predict_single(tspan,X,P,G,u,d,p)
% Fixed step Runge-Kutta 4 integrator
% Syntax: [Xk_RK4,Pk_RK4]=RK4_EKF_predict_single(tspan,X,P,G,u,d,p)
%         tspan: sampling time
%         X: state
%         P: covariance of the state
%         G: matrix for covariance
%         u: input
%         d: disturbance(meal) 
%         p: parameters of MVP model
M = 4; % RK4 steps per interval
DT = tspan/M;
Xk_RK4=X;
Pk_RK4=P;
for j=1:M
    [X_RK1,P_RK1] = CDEKFPredict_single(Xk_RK4, Pk_RK4, u, d, G,p);
    [X_RK2,P_RK2] = CDEKFPredict_single(Xk_RK4 + DT/2 * X_RK1, Pk_RK4 + DT/2 * P_RK1, u, d, G,p);
    [X_RK3,P_RK3] = CDEKFPredict_single(Xk_RK4 + DT/2 * X_RK2, Pk_RK4 + DT/2 * P_RK2, u, d, G,p);
    [X_RK4,P_RK4] = CDEKFPredict_single(Xk_RK4 + DT * X_RK3, Pk_RK4 + DT * P_RK3, u, d, G,p);
    Xk_RK4=Xk_RK4+DT/6*(X_RK1 +2*X_RK2 +2*X_RK3 +X_RK4);
    Pk_RK4=Pk_RK4+DT/6*(P_RK1 +2*P_RK2 +2*P_RK3 +P_RK4);
end
end