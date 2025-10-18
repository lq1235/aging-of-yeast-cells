clc;
clear;
close all
a1 = 0.1:0.4:18.0;
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
plot(a1, EPR, 'ks-', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2);
hold on;
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} EPR');
xlim([0,18])
ylim([0.179787894134981,0.367205853504862])
set(gca,'xtick',0:2:18)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);

h2=figure(2)
plot(a1, Flux, 'ko-', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2);
hold on;
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} Flux');
xlim([0,18])
ylim([0.06846685445993,0.11846685445993])
set(gca,'xtick',0:2:18)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);

h3=figure(3)
plot(a1,HDR,"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2);
hold on
xlabel('\fontsize{27} D');
ylabel('\fontsize{27} HDR');
xlim([0,18])
ylim([0.15,0.35])
set(gca,'xtick',0:2:18)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);

