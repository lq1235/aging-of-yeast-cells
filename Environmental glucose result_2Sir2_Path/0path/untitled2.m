clc;
clear;
Xm=2.5;
dif =[0.05;0.05];
ms = 200;
xmin = [0;0];
xmax = [Xm;Xm];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

px = load('pp_D1=1.10000.txt');
p = reshape(px(:,3),ms,ms);
sum(sum(p))
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;             
PP = eq(Pi,0)+Pi;    
P_eps=min(min(PP));   
P = P_eps*eq(Pi,0)+Pi;
U = -log(P);
U_max =3;
U=U.*(U<U_max)+U_max.*(U>U_max);
surf(X,Y,U)
hold on
shading interp;
colormap([jet(256)]);
view([-44.54,13.855]);
a=load('long-to-mode1.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.1,'m','LineWidth',3,'markersize',10)
hold on

a=load('mode1-to-long.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.1,'w','LineWidth',3,'markersize',10)
hold on

a=load('long-to-mode2.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.1,'m','LineWidth',3,'markersize',10)
hold on

a=load('mode2-to-long.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.1,'w','LineWidth',3,'markersize',10)
hold on

xlabel('\fontsize{25} Sir2');
ylabel('\fontsize{25} HAP');
zlabel('\fontsize{25} U')
axis([0 2,0 2])

set(gca,'xtick',0:0.5:3)
set(gca,'ytick',0:0.5:3)
set(gca,'ztick',0:2:4)
set(gca,'LineWidth',1,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
view([44.585472390627551,71.785433616115554])
