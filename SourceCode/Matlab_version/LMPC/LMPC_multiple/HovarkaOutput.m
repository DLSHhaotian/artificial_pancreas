function z = HovarkaOutput(x)
% hovorkaModel z = h(x) for Hovorka diabetic model
% Syntax: z = HovarkaOutput(x)
%         x: state
VG_BW=0.16;
BW=70;%[kg]
VG=VG_BW*BW;
z=x(5)/VG*18;%[mmol/L]-->[mg/dL]
end