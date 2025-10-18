clc
clear
a=load('BH.txt');
%MFPT
% h=figure(1)
% plot(a(:,1), a(:,2)-a(:,3), '^--', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
plot(a(:,1), a(:,2), 'p-', 'Color', 'm', 'LineWidth', 1, 'MarkerFaceColor', 'm', 'MarkerSize', 10)

% plot(a1, EPR, 'bs--', 'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 10);
% hold on
% plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
% hold on

% axis([-0.5 3.5 ,-1*10^(-3) 13*10^(-3)])
set(gca,'LineWidth',1.2,'Fontsize',24)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,-0.0035 0.04])
set(gca,'XTickLabelRotation',0);
ax = gca();
% ax.YRuler.Exponent = 0;
% ax.XRuler.Exponent = -2;
% set(gca,'xtick',0.056:0.0005:0.059)
xlabel("\fontsize{27} \alpha_{S0}");
ylabel("\fontsize{27}\Delta U");
xlim([0,1.2])
ylim([-0.5,5])
set(gca,'xtick',0:0.2:1.2)
set(gca,'ytick',0:1:5)
% ylim([-3.401295814688718,23.21870418531129])





% xlim([0.055173484275628,0.061162984275628])
% plot([0.7 0.7],[-0.5*10^(5) 9*10^(5)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-0.5*10^(5) 9*10^(5)],'b--','LineWidth',1)
% axis([-0.2 3.2,-0.5*10^(5) 9*10^(5)])



% %1/MFPT
% h=figure(2);
% plot(a(:,1),1./a(:,2),'ko-','LineWidth',1,'markersize',10)
% hold on
% % plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
% hold on

% set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
% set(gca,'LineWidth',1.2,'Fontsize',27)
% set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0 0.00395,-0.00005 0.001000001])
% set(gca,'xtick',0:0.0005:0.0040)
% xlabel("\fontsize{27} p");
% ylabel("\fontsize{27}\nu");
% xlim([0.001 0.0025])
% plot([0.7 0.7],[-1*10^(-4) 1.2*10^(-3)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-1*10^(-4) 1.2*10^(-3)],'b--','LineWidth',1)
% axis([-0.2 3.2,-1*10^(-4) 1.2*10^(-3)])
