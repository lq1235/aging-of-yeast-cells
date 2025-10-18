clc
clear
close all
h=figure(1);
a=load('data.txt');
plot(a(:,1),a(:,2),'b-', 'LineWidth', 1.5)
hold on
plot(a(:,1),a(:,3),'b-', 'LineWidth', 1.5)
hold on

b=xlsread('极限环数据.xlsx');
plot(b(2,:),b(1,:),'r-', 'LineWidth', 1.5)
hold on
plot(b(5,:),b(4,:),'black-', 'LineWidth', 1.5)
hold on

xlabel('\alpha_{H0}')
ylabel('S')
xlim([0 30])
ylim([0 400])

set(gca,'xtick',0:5:30)
set(gca,'ytick',0:100:400)
xtickangle(0);

set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
ax = gca;

