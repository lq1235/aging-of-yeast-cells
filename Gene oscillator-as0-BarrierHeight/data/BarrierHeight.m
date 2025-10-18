clc;
clear;
close all;
dif = [0.5; 0.5];
ms = 200;
xmin = [0;0];
xmax = [500;500];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);

[X,Y] = meshgrid(x,y);
b=load('UnstableNode.txt');    
a1 =0.1:0.1:1;

num = length(a1);
file_path = 'pp_D1=%0.5f.txt';
for j = 1:num
h=figure(j);
sample=sprintf(file_path,a1(j));
j
x1=b(j,1);
y1=b(j,2);
px = load(sample);
p = reshape(px(:,3),ms,ms);
sum(sum(p))
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;      
PP = eq(Pi,0)+Pi;     
P_eps=min(min(PP));  
P = P_eps*eq(Pi,0)+Pi;
P1=P/max(max(P));
U = -log(P1);
surf(X,Y,U)
shading interp;
colormap([jet(256)]);
hold on;
U_max=interp2(X,Y,U,x1,y1);
U_min=min(min(U))
BH=U_max-U_min
BH1(j,:)=[a1(j),BH];
set(gcf,'Position',[100 100 600 600]);
view([86.655,48.194]);
xlabel('\fontsize{25} Sir2');
ylabel('\fontsize{25} HAP');
zlabel('\fontsize{25} U');

axis([0 400,0 400])
set(gca,'xtick',0:100:500)
set(gca,'ytick',0:100:500)
% set(gca,'ztick',-10:5:18)
set(gca,'LineWidth',1.0,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [1 0.03])
set(gca,'XTickLabelRotation',0);
set(gca,'YTickLabelRotation',0);
end

