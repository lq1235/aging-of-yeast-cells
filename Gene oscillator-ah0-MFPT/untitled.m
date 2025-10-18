clc
clear
% a=load('MFPT.xlsx');
a=xlsread('MFPT.xlsx');
%MFPT
h=figure(1)

plot(a(:,1),a(:,2),'^-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
hold on
set(gca,'xtick',0:2:18)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
xlabel('\alpha_{H0}');
ylabel("\fontsize{27}T_{MFPT}");
ax = gca();
ax.YRuler.Exponent = 4;
set(gca,'XTickLabelRotation',0);
xlim([0 16])
ylim([-1914.943090185346,16235.05690981466])

%1/MFPT
h=figure(2);
plot(a(:,1),1./a(:,2),'o-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'XTickLabelRotation',0);
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0 0.00395,-0.00005 0.001000001])
set(gca,'xtick',0:2:18)
set(gca,'ytick',0.00:0.001:0.1)
xlabel('\alpha_{H0}');
ylabel("\fontsize{27}\nu");
xlim([0 16])
ylim([-0.000331834059515,0.004992165940485])
% plot([0.7 0.7],[-1*10^(-4) 1.2*10^(-3)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-1*10^(-4) 1.2*10^(-3)],'b--','LineWidth',1)
% axis([-0.2 3.2,-1*10^(-4) 1.2*10^(-3)])
