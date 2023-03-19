function xk = EulerMaruyamaSDE_MVP(tspan,dt,x,u,d,dw,p_sigma,p)
% Explicit Euler-Maruyama method with a fixed step size
% Syntax: xk = EulerMaruyamaSDE_MVP(tspan,dt,x,u,d,dw,p_sigma,p)
%         tspan: sampling time
%         dt: step size
%         x: state
%         u: input
%         d: disturbance(meal) 
%         dw: process noise
%         p_sigma: gain of process noise
%         p: parameters of MVP model
M = tspan/dt;
for i = 1:M
    xk = x;
    x = xk + MVPModel(tspan,x,u,d,p)*dt;
    x(4) = x(4) + p_sigma.*dw;
end
xk = x;
end