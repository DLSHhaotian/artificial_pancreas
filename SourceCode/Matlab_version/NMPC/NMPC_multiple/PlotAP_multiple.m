function PlotAP_multiple(x_arr,z_arr,y_arr,u_arr,d_arr,HR_arr,Z_bar,t_cost)
t_hour=0:minutes(5):minutes(24*60);
fs = 9; % Font size
f=gobjects(0);
f(end+1)=figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);

%figure;
width=900;
height=1000;
left=0;
bottem=0;
set(gcf,'position',[left,bottem,width,height])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(22,1,1:4)
fill([minutes(0),minutes(1440),minutes(1440),minutes(0)],[70,70,140,140],'g')
hold on;
plot(t_hour,z_arr,'LineWidth',3.0)
title('Blood Glucose')
%legend("Normal range","CGM blood glucose","Location",'best');
ylabel('G [mg/dL]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%meal
subplot(22,1,6:7)
stem(t_hour(1:end-1),d_arr*5,'LineWidth',3.0,'MarkerSize',4)
title('Size of the meal')
%legend("Disturbance d","Location", 'best');
%xlabel('Time [h]');
ylabel('D [g]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%xlh = xlabel('Time[h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HR
subplot(22,1,9:10)
stairs(t_hour(1:end-1),HR_arr,'LineWidth',3.0,'MarkerSize',4)
title('Excercise')
%legend("Disturbance d","Location", 'best');
%xlabel('Time [h]');
ylabel('HR [bpm]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%xlh = xlabel('Time[h]');
%xlh.Position(2) = xlh.Position(2) - 5.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(22,1,12:13)
stairs(t_hour(1:end-1)',u_arr(1,:)',':','LineWidth',3.0)
title('Basal insulin')
ylabel('Uba [mU/min]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(22,1,15:16)
stem(t_hour(1:end-1)',u_arr(2,:)',':','LineWidth',3.0,'MarkerSize',4)
title('Bolus insulin')
ylabel('Ubo [mU/min]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%uglu
subplot(22,1,18:19)
stem(t_hour(1:end-1),u_arr(3,:)',':','LineWidth',3.0,'MarkerSize',4)
title('Glucagon')
ylabel('Uglu [{\mu}g]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tcost
subplot(22,1,21:22)
stairs(t_hour(1:end-1)',t_cost(1,:)',':','LineWidth',3.0,'MarkerSize',4)
title('Calculation time')
xlabel('Time [h]');
ylabel('Cost [s]');
xtickformat('hh:mm')
set(gca,'fontsize',fs)
set(gca,'XTick',[t_hour(1):minutes(6*60):t_hour(end)]) %

end