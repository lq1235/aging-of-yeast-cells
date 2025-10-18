clc
clear
a=load('MFPT_h=0.002_r=0.01.txt');

h=figure(1)
plot(a(:,1),a(:,2),'^-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 6.2)
hold on
xlim([0 5])
ylim([21.141003231775237,118.6195211265201])
set(gca,'XTickLabelRotation',0);
set(gca,'xtick',0:1:5)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.1])
xlabel("\fontsize{27} Glucose%");
ylabel("\fontsize{27}T_{MFPT}");

h=figure(2);
plot(a(:,1),1./a(:,2),'^-', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 6.2)
hold on
xlim([0 5])
ylim([0.005432963104116,0.039315756553876])
set(gca,'xtick',0:1:5)
set(gca,'XTickLabelRotation',0);
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.1])
xlabel("\fontsize{27} Glucose%");
ylabel("\fontsize{27}\nu");

