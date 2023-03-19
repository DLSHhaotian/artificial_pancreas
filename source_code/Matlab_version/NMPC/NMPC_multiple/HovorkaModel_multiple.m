function xdot = HovorkaModel_multiple(t,x,u,d,p,HR)
% hovorkaModel dx/dt = f(t,x,u,d,p) for Hovorka diabetic model
%
% This function implements a differential equation model for the
% Hovorka diabetic model.
%
% Syntax: xdot = HovorkaModel_single(t,x,u,d,p)
%         x: gradient of x
%         u: steady-state of input uss
%         y: steady-state of output yss
%         d: disturbance(meal) steady-state is 0
%         p: parameters of Hovorka model
%         HR: disurbance(Heart rate)

%%%%%model parameters%%%%%
%MwG: [g/mol] the molecular weight of glucose
%AG: the carbohydrate bioavailability parameter
%tau_D: [min] the glucose absorption rate(time constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VG: [L] the glucose distribution volume
%F01: [mmol/min] the glucose consumption of the central nervous system at high glucose concentration
%k12: [1/min] the transfer rate from the blood to the tissues
%EGP0 :[mmol/min] the amount of glucose produced by the liver at zero insulin level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tau_S: [min] the subcutaneous-to-intravenous absorption rate(time constant)
%VI: [L] the insulin distribution volume
%ke: [1/min] is the insulin elimination rate
%ka1: [1/min] constants
%ka2: [1/min] constants
%ka3: [1/min] constants
%SIT: [L/mU] insulin sensitivity
%SID: [L/mU] insulin sensitivity
%SIE: [L/mU] insulin sensitivity
%kb1: [L/mU/min] SIT*ka1
%kb2: [L/mU/min] SID*ka2
%kb3: [L/mU/min] SIE*ka3
%%%%%Input%%%%%%
%x=[x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11]
%u=[uba; ubo; dmeal]
%d=[dw]
%p=[MwG; AG; tau_D; VG; F01; k12; EGP0; tau_S; VI; ke; ka1; ka2; ka3; SIT; SID; SIE; tau_GI]

p1=p(1:3);%Food Absorption
p2=p(4:7);%Glucose
p3=p(8:end);%Insulin
uba=u(1);%basal insulin
ubo=u(2);%bolus insulin
uglu=u(3);%glucagon
dmeal=d;%meal intake

%HR = parm.HR; % Heart rate [bpm]
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
xdot(12)= uglu-x(12)/tauGlu;
xdot(13)= (x(12)-x(13))/tauGlu;

% Exercise subsystem
dHR = HR-HRbase;
xdot(14) = 1/tHR*(dHR-x(14));
fE1 = ((x(14)/(a*HRbase))^n)/(1+(x(14)/(a*HRbase))^n);
xdot(15) = 1/tex * (c1*fE1+c2-x(15));
xdot(16) = -(fE1/tin+1/x(15))*x(16) + (fE1*x(15))/(c1+c2);
end