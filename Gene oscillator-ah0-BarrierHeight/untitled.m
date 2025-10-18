clc
clear
a=load('BH.txt');
plot(a(:,1), a(:,2), 'p-', 'Color', 'm', 'LineWidth', 1, 'MarkerFaceColor', 'm', 'MarkerSize', 10)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);
ax = gca();
set(gca,'xtick',0:2:16)
xlabel("\fontsize{27} \alpha_{H0}");
ylabel("\fontsize{27}\Delta U");
xlim([0,16])
ylim([-0.517838455596066,5])
set(gca,'xtick',0:2:16)
set(gca,'ytick',0:1:8)


