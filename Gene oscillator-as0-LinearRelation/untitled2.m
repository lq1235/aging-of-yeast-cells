clc
clear
a=load('barrier_height.txt');


% plot(a(:,2)-a(:,3),log(b(:,2)),'x-k','LineWidth',1,'Markersize',10);

plot(a(:,2),log(a(:,3)),'o-', 'Color', 'black', 'LineWidth', 1, 'MarkerFaceColor', 'black', 'MarkerSize', 10);


% axis([0 15,3 16])
xlabel('\fontsize{27}\DeltaU')
ylabel('\fontsize{27}logT_{MFPT}')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'xtick',0:1:6)
set(gca,'ytick',0:1:10)
xlim([-0.422241708649475,4.901758291350524])
ylim([4.623618759132745,10.57588845731484])