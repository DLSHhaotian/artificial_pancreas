function x = HovorkaModely2x_single(u,y,d,p)
%Use after calculating the steady-state output and input to calculate the steady-state xss

% Function for Hovorka model 
% Syntax: x = HovorkaModely2x_single(u,y,d,p)
%         x: steady state xss
%         u: steady-state of input uss
%         y: steady-state of output yss
%         d: disturbance(meal) steady-state is 0
%         p: parameters of MVP model


p1=p(1:3);%Food Absorption
p2=p(4:7);%Glucose
p3=p(8:end);%Insulin
uba=u(1);%basal insulin
ubo=u(2);%bolus insulin
dmeal=d;%meal intake
%p=[MwG; AG; tau_D;   VG; F01; k12; EGP0;   tau_S; VI; ke; ka1; ka2; ka3; SIT; SID; SIE;tau_GI]
x=zeros(11,1);
x(1)=0;
x(2)=0;
x(3)=uba*p3(1);
x(4)=uba*p3(1);
x(5)=p2(1)*y;
x(6)=(p2(1)*y*p3(7)*p3(4)*uba)/(p3(4)*p3(2)*p3(3)*(p2(3)+p3(8)*p3(5)*uba/(p3(5)*p3(2)*p3(3))));
x(7)=uba/(p3(2)*p3(3));
x(8)=(p3(7)*p3(4)*uba/(p3(2)*p3(3)))/p3(4);
x(9)=(p3(8)*p3(5)*uba/(p3(2)*p3(3)))/p3(5);
x(10)=(p3(9)*p3(6)*uba/(p3(2)*p3(3)))/p3(6);
x(11)=y;
end