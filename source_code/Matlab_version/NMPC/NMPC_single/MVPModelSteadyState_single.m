function xdot = MVPModelSteadyState_single(u,y,d,p)
% Use to calculate the steady-state output and input(using fsolve function)

% Function for MVP model 
% Syntax: xdot = MVPModelSteadyState_single(u,y,d,p)
%         x: gradient of x
%         u: steady-state of input uss
%         y: steady-state of output yss
%         d: disturbance(meal) steady-state is 0
%         p: parameters of MVP model


%MVP p=[k1; CI; SI; EGP; km; VG; tau_GI; rou4]
uba=u(1);%basal insulin
ubo=u(2);%bolus insulin
dmeal=d;%meal intake
x=zeros(7,1);
x(1)=uba/p(2);
x(2)=uba/p(2);
x(3)=uba*p(3)/p(2);
x(4)=y;
x(5)=0;
x(6)=0;
x(7)=y;

% Differential equations
xdot = zeros(7,1);
xdot(1,1)=p(1)*((uba+ubo)/p(2)-x(1));
xdot(2,1)=p(1)*(x(1)-x(2));
xdot(3,1)=-p(1)*x(3)+p(1)*p(3)*x(2);
xdot(4,1)=-x(3)*x(4)+p(4)+(x(6)*p(5)/p(6));
xdot(5,1)=dmeal-x(5)*p(5);
xdot(6,1)=(x(5)-x(6))*p(5);
xdot(7,1)=(x(4)-x(7))/p(7);

end