function [hdot,z] = MVPOutputJacobian_single(x)
%  [hdot,z] = jacobian(h(x)) for MVP model
% Syntax: [hdot,z] = MVPOutputJacobian_single(x)
%         hdot: jacobian matrix of h(x)
%         z: output
%         x: current state

hdot=zeros(1,7);
hdot(4)=1;
z=x(4);
end