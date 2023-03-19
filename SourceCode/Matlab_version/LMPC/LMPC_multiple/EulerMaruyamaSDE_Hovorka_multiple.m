function xk = EulerMaruyamaSDE_Hovorka_multiple(tspan,dt,x,u,d,dw,p_sigma,p,HR)
% Explicit Euler-Maruyama method with a fixed step size
% Syntax: xk = EulerMaruyamaSDE_Hovorka_multiple(tspan,dt,x,u,d,dw,p_sigma,p,HR)
%         tspan: sampling time
%         dt: step size
%         x: state
%         u: input
%         d: disturbance(meal) 
%         dw: process noise
%         p_sigma: gain of process noise
%         p: parameters of Hovorka model
%         HR: heart rate

M = tspan/dt;
for i = 1:M
    xk = x;
    x = xk + HovorkaModel_multiple(dt,xk,u,d,p,HR)*dt + p_sigma.*dw;
end
xk = x;
end