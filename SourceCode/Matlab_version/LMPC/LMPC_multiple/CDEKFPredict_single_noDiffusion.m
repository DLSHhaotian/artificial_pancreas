function [Xk_pred,Pk_pred] = CDEKFPredict_single_noDiffusion(X,P,u,d,G,p)
% Syntax: [Xk_pred,Pk_pred] = CDEKFPredict_single_noDiffusion(X,P,u,d,G,p)
%         X: state
%         P: covariance of the state
%         u: input
%         d: disturbance(meal) 
%         G: matrix for covariance
%         p: parameters of MVP model

% --------------------------------------------------------------
% Prediction of state
% --------------------------------------------------------------
[fdot,Xk_pred] = MVPStateJacobian_single_noDiffusion(X,u,d,p);
% --------------------------------------------------------------
% Prediction of covariance
% --------------------------------------------------------------
Pk_pred=fdot*P+P*fdot'+G*G';
end