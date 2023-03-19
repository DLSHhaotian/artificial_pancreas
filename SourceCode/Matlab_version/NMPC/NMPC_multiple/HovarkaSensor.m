function y = HovarkaSensor(x,noise)
% hovorkaModel y = g(x) for Hovorka diabetic model
% Syntax: y = HovarkaSensor(x,noise)
%         x: state
%         noise: sensor noise
%x(11): CGM[mmol/L]
%y: CGM[mg/dL]
y=x(11)*18+noise;
end

