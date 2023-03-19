function xdot = HovorkaModelSteadyState_multiple(u,y,d,p,HR)
%Use to calculate the steady-state output and input(using fsolve function)
% Syntax: xdot = HovorkaModelSteadyState_multiple(u,y,d,p,HR)
%         x: gradient of x
%         u: steady-state of input uss
%         y: steady-state of output yss
%         d: disturbance(meal) steady-state is 0
%         p: parameters of Hovorka model
%         HR: disurbance(Heart rate)


a = 0.77; %Exercise parameter [unitless]
tHR = 5; %Exercise parameter [min]
tin = 1; %Exercise parameter [min]
n = 3; %Exercise parameter [unitless]
tex = 200; %Exercise parameter [min]
c1 = 500; %Exercise parameter [min]
c2 = 100; %Exercise paramtere [min]
beta = 0.78; %Exercise-induced insulin-independent glucose uptake rate [mmol/min]
alpha = 10;%1.79; % Exercise-induced insulin action [unitless]
HRbase = 70; %Basal heart rate [BPM]

tauGlu = 19; % Glucagon time constant [min]
KGlu = 0.075; % Glucagon gain

p1=p(1:3);%Food Absorption
p2=p(4:7);%Glucose
p3=p(8:end);%Insulin
uba=u(1);%basal insulin
ubo=u(2);%bolus insulin
uglu=u(3);%glucagon
dmeal=d;%meal intake
%p=[MwG; AG; tau_D;   VG; F01; k12; EGP0;   tau_S; VI; ke; ka1; ka2; ka3; SIT; SID; SIE;tau_GI]
x=zeros(16,1);
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
x(12)=0;
x(13)=0;
x(14)=0;
x(15)=c2;
x(16)=0;
G=x(5)/p2(1);%glucose concentration
%The glucose consumption of the central nervous system
if G>=4.5
    F01_c=p2(2);
else
    F01_c=p2(2)*G/4.5;
end
%The production of glucose in the kidneys 
if G>=9
    FR=0.003*(G-9)*p2(1);
else
    FR=0;
end


%p=[MwG; AG; tau_D;   VG; F01; k12; EGP0;   tau_S; VI; ke; ka1; ka2; ka3; SIT; SID; SIE;tau_GI]
% Differential equations
xdot = zeros(16,1);
%meal
Dt=1000*dmeal/p1(1);
xdot(1,1) = p1(2)*Dt-x(1)/p1(3); 
xdot(2,1) = x(1)/p1(3)-x(2)/p1(3);
UG=x(2)/p1(3);
%glucose
xdot(5,1) = UG-F01_c-FR-x(8)*x(5)+p2(3)*x(6)+p2(4)*(1-x(10))+ KGlu*p2(1)*x(13)/18.01577 - (alpha*x(16)^2)*x(8)*x(5); 
xdot(6,1) = x(8)*x(5)-(p2(3)+x(9))*x(6)+ (alpha*x(16)^2)*x(8)*x(5) - x(9)*x(6)*(alpha*x(16)^2) - beta*x(14)/HRbase;
%insulin
UI=x(4)/p3(1);
xdot(3,1) = uba+ubo-x(3)/p3(1); 
xdot(4,1) = x(3)/p3(1)-x(4)/p3(1);
xdot(7,1) = UI/p3(2)-p3(3)*x(7); 
xdot(8,1) = -p3(4)*x(8)+p3(7)*p3(4)*x(7);
xdot(9,1) = -p3(5)*x(9)+p3(8)*p3(5)*x(7);
xdot(10,1) = -p3(6)*x(10)+p3(9)*p3(6)*x(7);
%CGM
xdot(11,1) = G/p3(10)-x(11)/p3(10);
% Subcutaneous glucagon
xdot(12,1)= uglu-x(12)/tauGlu;
xdot(13,1)= (x(12)-x(13))/tauGlu;

% Exercise subsystem
dHR = HR-HRbase;
xdot(14,1) = 1/tHR*(dHR-x(14));
fE1 = ((x(14)/(a*HRbase))^n)/(1+(x(14)/(a*HRbase))^n);
xdot(15,1) = 1/tex * (c1*fE1+c2-x(15));
xdot(16,1) = -(fE1/tin+1/x(15))*x(16) + (fE1*x(15))/(c1+c2);
end