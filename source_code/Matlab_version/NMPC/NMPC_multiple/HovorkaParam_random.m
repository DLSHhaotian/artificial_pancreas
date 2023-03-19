function p_Hovorka = HovorkaParam_random(seed)
% Syntax: p_Hovorka = HovorkaParam_random(seed)
%         seed: seed for random number 

% Control random number 
rng(seed)
% Number of uniformly distributed parameters
n_uniform = 2; 
% Number of normally distributed parameters
n_normal = 14; 

% Uniform distribution
arr_uniform = rand(n_uniform,1); 
% Normal distribution
arr_normal = randn(n_normal,1);  
% Normal distribution
while(sum(arr_normal>1)~=0 || sum(arr_normal<-1)~=0)
    arr_normal = randn(n_normal,1);  
end

% Generate parameters
BW = 80+30*(arr_uniform(1)-0.5);
EGP0 = (0.0161+1.6e-3*sqrt(6)*arr_normal(1))*BW;
F01 = (0.0097+0.9e-3*sqrt(6)*arr_normal(2))*BW;
Ag = 0.7+arr_uniform(2)*(1.2 - 0.7);
k12 = 0.0649+sqrt(6)*0.0115*arr_normal(3);
ka1 = 0.0055+sqrt(6)*0.0023*arr_normal(4);
ka2 = 0.0683+sqrt(6)*0.0207*arr_normal(5);
ka3 = 0.0304+sqrt(6)*0.0096*arr_normal(6);
SIT = (51.2+13.1*sqrt(6)*arr_normal(7))*1e-4;
SID = (8.2+3.2*sqrt(6)*arr_normal(8))*1e-4;
SIE = (520+125*sqrt(6)*arr_normal(9))*1e-4;
ke = 0.14+0.035*arr_normal(10);
Vi = (0.12+0.012*arr_normal(11))*BW;
Vg = exp(log(0.15)+0.23*arr_normal(12))*BW;
tmaxI = 1/(0.018+0.0045*arr_normal(13));
tmaxG = 1/(exp(-3.689+0.25*arr_normal(14)));

%Hovorka p=[MwG; AG; tau_D; VG; F01; k12; EGP0; tau_S; VI; ke; ka1; ka2; ka3; SIT; SID; SIE; tau_GI]
%tau_D = tau_G;tau_S = tau_I;
MwG=180.1577;%[g/mol]
tau_GI=6.7;
p_Hovorka = [MwG; Ag; tmaxG; Vg; F01; k12; EGP0; tmaxI; Vi; ke; ka1; ka2; ka3; SIT; SID; SIE; tau_GI];
end