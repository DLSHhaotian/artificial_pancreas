function p_MVP = ParamEstimateMLE(t_arr,u_arr,d_arr,y_arr,dW,T,p_MVP,xs_MVP)
% Syntax: p_MVP = ParamEstimateMLE(t_arr,u_arr,d_arr,y_arr,dW,T,p_MVP,xs_MVP)
%         t_arr: time of data set
%         u_arr: input of data set
%         d_arr: disturbance of data set
%         y_arr: output of data set
%         dW: process noise of data set
%         T: sampling time
%         p_MVP: parameters of MVP model
%         xs_MVP: steady state of MVP model

uba_data=u_arr(1,:);
ubo_data=u_arr(2,:);
d_data=d_arr;
y_data=y_arr;
t_data=t_arr';

p0=p_MVP;
x0=xs_MVP;

n_MVP = size(p0,1);
lb = zeros(n_MVP,1);
options = optimset('display','final',...
'MaxIter', 5000,...
'Algorithm','sqp');

xsol = fmincon(@(x) objfunMLE(x,t_data, y_data, [uba_data;ubo_data],d_data*1000,dW./T, x0),p0,[],[],[],[],lb,[],[],options);
p_MVP = xsol;
end