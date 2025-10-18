clc
clear
% a=load('MFPT.xlsx');
a=xlsread('MFPT.xlsx');
%MFPT
h=figure(1)
plot(a(:,1),a(:,2),'^-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
hold on
% plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
% hold on
set(gca,'xtick',0:0.2:1.2)
set(gca,'LineWidth',1.2,'Fontsize',24)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
xlabel('\alpha_{S0}');
ylabel("\fontsize{27}T_{MFPT}");
ax = gca();
ax.YRuler.Exponent = 4;
set(gca,'XTickLabelRotation',0);
xlim([0 1.2])
ylim([-2463.205867392261,19536.79413260775])
% axis([-0.2 3.2,-0.5*10^(5) 9*10^(5)])
% xlim([0.001 0.0025])


%1/MFPT
h=figure(2);
plot(a(:,1),1./a(:,2),'o-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
hold on
% plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
% hold on

set(gca,'LineWidth',1.2,'Fontsize',24)
set(gca,'XTickLabelRotation',0);
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0 0.00395,-0.00005 0.001000001])
set(gca,'xtick',0:0.2:1.2)
xlabel('\alpha_{S0}');
ylabel("\fontsize{27}\nu");
xlim([0 1.2])
ylim([-0.001034232870032,0.008645767129968])
% plot([0.7 0.7],[-1*10^(-4) 1.2*10^(-3)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-1*10^(-4) 1.2*10^(-3)],'b--','LineWidth',1)
% axis([-0.2 3.2,-1*10^(-4) 1.2*10^(-3)])
