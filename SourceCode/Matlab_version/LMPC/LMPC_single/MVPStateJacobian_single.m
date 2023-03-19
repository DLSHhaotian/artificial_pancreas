function [fdot,xdot] = MVPStateJacobian_single(x,u,d,p)
%  [fdot,xdot] = jacobian(f(x)) for MVP model

% Syntax: [fdot,xdot] = MVPStateJacobian_single(x,u,d,p) 
%         fdot: jacobian matrix of f(x)
%         xdot: gradient of x
%         x: current state
%         u: current input[uba;ubo;dmeal]
%         d: current process noise
%         p: parameters of MVP model

fdot=zeros(7,7);
fdot(1,1)=-p(1);
fdot(2,1)=p(1);
fdot(2,2)=-p(1);
fdot(3,2)=p(3)*p(1);
fdot(3,3)=-p(1);
fdot(4,3)=-x(4);
fdot(4,4)=-x(3);
fdot(4,6)=p(5)/p(6);
fdot(5,5)=-p(5);
fdot(6,5)=p(5);
fdot(6,6)=-p(5);
fdot(7,4)=1/p(7);
fdot(7,7)=-1/p(7);
%MVP p=[k1; CI; SI; EGP; km; VG; tau_GI; rou4]
uba=u(1);%basal insulin
ubo=u(2);%bolus insulin
dmeal=u(3);%meal intake
dw=d;%disturbance


% Differential equations
xdot = zeros(7,1);
%xdot(1,1)=p(1)*((uba+ubo)/p(2)-x(1));
xdot(1,1)=p(1)*((uba+ubo)/1.5-x(1));
xdot(2,1)=p(1)*(x(1)-x(2));
xdot(3,1)=-p(1)*x(3)+p(1)*p(3)*x(2);
xdot(4,1)=-x(3)*x(4)+p(4)+(x(6)*p(5)/p(6))+p(8)*dw;
xdot(5,1)=dmeal-x(5)*p(5);
xdot(6,1)=(x(5)-x(6))*p(5);
%xdot(7,1)=(x(4)-x(7))/p(7);
xdot(7,1)=(x(4)-x(7))/6.7;

end