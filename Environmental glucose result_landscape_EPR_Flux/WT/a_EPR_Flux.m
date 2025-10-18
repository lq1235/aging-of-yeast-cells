clc;
clear;
close all
a1 = 0.01:0.03:0.5;
dif =[0.05;0.05];
ms = 200;
xmin = [0;0];
xmax = [2.5;2.5];
num = length(a1);
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);
file_path = 'pp_D1=%0.5f.txt';
for j = 1:num
    sample=sprintf(file_path,a1(j));
    px = load(sample);
    p = reshape(px(:,3),ms,ms);
    FPx = reshape(px(:,4),ms,ms);
    FPy = reshape(px(:,5),ms,ms);
    z = trapz(y,trapz(x,p));
    Pi = p/z;
    PP = eq(Pi,0)+Pi;
    P_eps=min(min(PP));
    P = P_eps*eq(Pi,0)+Pi;
    U = -log(P);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    [GPx,GPy] = gradient(P,dx,dy);
    Jx = FPx.*P - dif(1)*GPx ;
    Jy = FPy.*P - dif(2)*GPy ;
    JJP=(Jx.^2)./(dif(1)*P)+(Jy.^2)./(dif(2)*P);
    EPR(j)=trapz(y,trapz(x,JJP));
    FJ=(FPx.*Jx+FPy.*Jy)/dif(1);
    HDR(j)=trapz(y,trapz(x,FJ));
    JJ=sqrt(Jx.^2+Jy.^2);
    Flux(j) = trapz(y,trapz(x,JJ));
end

h1=figure(1)
plot(a1, EPR, 'bs-', 'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 6.2);
hold on;
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} EPR');
xlim([0,0.5])
ylim([0.226821442249296,0.27801011066377])
set(gca,'xtick',0:0.1:0.5)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
ax = gca;
ax.YAxis(1).Exponent = -1;

h2=figure(2)
plot(a1, Flux, 'ro-', 'LineWidth', 1, 'MarkerFaceColor', 'red', 'MarkerSize', 6.2);
hold on;
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} Flux');
xlim([0,0.5])
ylim([0.093435565442284,0.102235565442284])
set(gca,'xtick',0:0.1:0.5)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
ax = gca;
ax.YAxis(1).Exponent = -2;
print(h2, '-r600', '-dpdf', ['Flux', num2str(1),'.pdf']);

h3=figure(3)
plot(a1,HDR,"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2);
hold on
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} HDR');
xlim([0,0.5])
ylim([0.200725990296392,0.244725990296392])
ax = gca;
ax.YAxis(1).Exponent = -1;
set(gca,'xtick',0:0.1:0.5)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])