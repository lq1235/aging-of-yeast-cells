clc
clear
close all
a=load('势垒高度.txt');
plot(a(:,1), a(:,3)-a(:,2), 'bo-', 'Color', 'black', 'LineWidth', 1, 'MarkerFaceColor', 'black', 'MarkerSize', 6.2)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);
ax = gca();
ax.YRuler.Exponent = 0;
set(gca,'xtick',0:1:5)
xlabel("\fontsize{27} Glucose%");
ylabel("\fontsize{27}\Delta U");
xlim([0,5])
ylim([-0.2,1.5])

