clc
clear
a=xlsread('data.xlsx');

plot(a(:,2),a(:,3),"black-o",'LineWidth', 1, 'MarkerFaceColor', 'black', 'MarkerSize', 10)
xlabel('\fontsize{24} Coherence');
ylabel('\fontsize{24} Flux');

% ax = gca();
% ax.YRuler.Exponent = -1;
% ax.XRuler.Exponent = -1;
% set(gca,'xtick',0:1:15)
% xlim([0,1.2])
% ylim([0.93260142627323,1.006174457294026])
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
% axis([0 0.00395,0 0.00025])
set(gca,'XTickLabelRotation',0);

