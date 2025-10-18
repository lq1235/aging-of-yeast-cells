
clc
clear
figure(1)
a=xlsread('去S较小的态.xlsx');
yyaxis left
h1=plot(a(1:2:end,1), a(1:2:end,3)-a(1:2:end,2),'ro-', 'LineWidth', 1, 'MarkerFaceColor', 'red', 'MarkerSize', 6.2);
hold on
xlim([-0.3   5]);
ylim([  0.38   1.21]);
ylabel('\fontsize{27} \Delta U');
b=xlsread('去S较大的态.xlsx');
yyaxis right
h2=plot(b(1:2:end,1), b(1:2:end,3)-b(1:2:end,2),'bo-', 'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 6.2);
xlim([-0.3    5]);
ylim([ -0.05    1.7]);
leg = legend([h1, h2], 'Mode 1', 'Mode 2', 'fontsize', 10);
set(leg, 'Position', [0.4, 0.8, 0.2, 0.1]);
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'ycolor','g');
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';
set(gca,'xtick',0:1:5)
set(gca,'XTickLabelRotation',0);
xlabel('\fontsize{27} Glucose%');
ylabel('\fontsize{27} \Delta U');
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])

