function PlotAP_single(x_arr,z_arr,y_arr,u_arr,d_arr,Z_bar,t_cost)
t_hour=0:minutes(5):minutes(24*60);
fs = 11; % Font size
f=gobjects(0);
f(end+1)=figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);

%figure;
width=1400;
height=1500;
left=0;
bottem=0;
set(gcf,'position',[left,bottem,width,height])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(19,1,1:4)
fill([minutes(0),minutes(1440),minutes(1440),minutes(0)],[70,70,140,140],'g')
hold on;
plot(t_hour,z_arr,'LineWidth',3.0)
title('Blood Glucose')
legend("Normal range","CGM blood glucose","Location",'best');
%xlh = xlabel('Time [h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;
ylabel('G [mg/dL]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(19,1,6:7)
stem(t_hour(1:end-1),d_arr*5,'LineWidth',3.0,'MarkerSize',2)
title('Size of the meal')
%legend("Disturbance d","Location", 'best');
%xlabel('Time [h]');
ylabel('D [g]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%xlh = xlabel('Time[h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(19,1,9:11)
%plot(t_hour(1:end-1),u_duba(1,:),t_hour(1:end-1),u_ubabar(1,:),t_hour(1:end-1),u_duba_ubo(1,:),t_hour(1:end-1),u_ubabar_ubo(1,:),'LineWidth',3.0)
stairs(t_hour(1:end-1)',u_arr(1,:)',':','LineWidth',3.0)
%hold on;
%stairs(t_hour(1:end-1)',u_duba_ubo(1,:)',':','LineWidth',3.0)
%legend("Matlab version","C++ version","Location",'best');

title('Basal insulin')
%xlabel('Time [h]');
ylabel('Uba [mU/min]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%xlh = xlabel('Time[h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(19,1,13:15)
%stairs(t_hour(1:end-1)',[u_duba(2,:)',u_ubabar(2,:)',u_ubabar_ubo(2,:)',u_duba_ubo(2,:)'],':','LineWidth',3.0)
stem(t_hour(1:end-1)',u_arr(2,:)',':','LineWidth',3.0,'MarkerSize',4)
%legend("Matlab version","C++ version","Location",'best');

title('Bolus insulin')
%xlabel('Time [h]');
ylabel('Ubo [mU/min]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%xlh = xlabel('Time[h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(19,1,17:19)
%stairs(t_hour(1:end-1)',[u_duba(2,:)',u_ubabar(2,:)',u_ubabar_ubo(2,:)',u_duba_ubo(2,:)'],':','LineWidth',3.0)
stairs(t_hour(1:end-1)',t_cost(1,:)',':','LineWidth',3.0,'MarkerSize',4)
%legend("Matlab version","C++ version","Location",'best');

title('Calculation time')
xlabel('Time [h]');
ylabel('Cost [s]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
end