clc
clear
close all
h=figure(1);
a=load('data.txt');
plot(a(:,1),a(:,2),'b-', 'LineWidth', 1.5)
hold on
plot(a(:,1),a(:,3),'b-', 'LineWidth', 1.5)
hold on

b=xlsread('data2.xlsx');
plot(b(:,1),b(:,2),'r-', 'LineWidth', 1.5)
hold on

c=load('data1.txt');
plot(c(:,1),c(:,2),'black-', 'LineWidth', 1.5)
hold on

xlabel('\alpha_{S0}')


ylabel('S')
xlim([0 5])
ylim([0 400])


set(gca,'xtick',0:1:5)
set(gca,'ytick',0:100:400)
xtickangle(0);

set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
ax = gca;
print(h, '-r600', '-dpdf', ['Phase', num2str(1),'.pdf']);
