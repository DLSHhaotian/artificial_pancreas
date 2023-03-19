function x = MVPModely2x_single(u,y,d,p)
% Use after calculating the steady-state output and input to calculate the steady-state xss

% Function for MVP model 
% Syntax: x = MVPModely2x_single(u,y,d,p)          
%         x: steady state xss
%         u: steady-state of input uss
%         y: steady-state of output yss
%         d: disturbance(meal) steady-state is 0
%         p: parameters of MVP model
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
end