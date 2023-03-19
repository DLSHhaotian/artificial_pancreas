function obj = objfunMLE(x,t_data, y_data, u_data, d_data, dW, x0)
% Syntax: obj = objfunMLE(x,t_data, y_data, u_data, d_data, dW, x0)
%         x: parameters of MVP model
%         t_data: time of data set
%         y_data: output of data set
%         u_data: input of data set
%         d_data: disturbance of data set
%         dW: process noise of data set
T=5;

%MVP p=[k1; CI; SI; EGP; km; VG; tau_G; tau_GI; rou4]
len=size(t_data,2);
%weight=(1./y_data).^2;
%weight = ones(1,len);

xk=zeros(1,len);
xk(1,1)=x0(4);
obj=(len+1)*len/2*log(2*3.14);
x_bar=x0;

G=zeros(7,1);
G(4)=x(8);
P_kf = G*G';
%P_kf=eye(7);
R_kf = 0.02^2;

for k=1:len-1
    %%%%update%%%%
    [C_kf,z_bar] = MVPOutputJacobian_single(x_bar);
    e_kf=y_data(:,k)-z_bar;
    R_kf=C_kf*P_kf*C_kf'+R_kf;
    K_kf=P_kf*C_kf'/R_kf;
    x_bar=x_bar+K_kf*e_kf;
    P_kf=P_kf-K_kf*R_kf*K_kf';
   
    %%%%predict%%%%
    [x_bar,P_kf]=EulerMaruyama_EKF_predict_single(T,x_bar,P_kf,G,[u_data(:,k);d_data(:,k)],dW(k),x);%MVP model meal[mg/min]
 
    obj=obj+0.5*(log(R_kf)+e_kf'/R_kf*e_kf);
end