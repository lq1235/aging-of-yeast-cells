clc;
clear;
dif = [0.05; 0.05];
ms = 200;
xmin = [0;0];
xmax = [2.5;2.5];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);
a1 =0.01:0.03:0.01;
num = length(a1);
file_path = 'pp_D1=%0.5f.txt';
for j = 1:num
    h=figure(j);
    sample=sprintf(file_path,a1(j));
    px = load(sample);
    p = reshape(px(:,3),ms,ms);
    sum(sum(p));
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
    hold on
    zlim([-5 15])
    view([-40.422194513715709,9.402001279491998]);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    [GPx,GPy] = gradient(P,dx,dy);
    Jx = FPx.*P - dif(1)*GPx ;
    Jy = FPy.*P - dif(2)*GPy ;
    E=Jy.^2+Jx.^2;
    JJx=Jx./(sqrt(E)+eps);
    JJy=Jy./(sqrt(E)+eps);
    Fx= dif(1)*GPx./P;
    Fy=dif(2)*GPy./P;
    F=Fx.^2+Fy.^2;
    FFx=Fx./(sqrt(F)+eps);
    FFy=Fy./(sqrt(F)+eps);
    mg = 5:13:200;
    ng = mg;
    % quiver(X(mg,ng),Y(mg,ng),JJx(mg,ng),JJy(mg,ng),0.6,'color','w','LineWidth',1);
    hold on
    % quiver(X(mg,ng),Y(mg,ng),FFx(mg,ng),FFy(mg,ng),0.5,'color','k');
    axis([0 2.5,0 2.5])
    xlabel('\fontsize{35} Sir2');
    ylabel('\fontsize{35} HAP');
    zlabel('\fontsize{35} U');
    set(gca,'xtick',-0:1:2.5)
    set(gca,'ytick',-0:1:2.5)
    set(gca,'LineWidth',1.0,'Fontsize',35)
    set(gca,'TickDir', 'out', 'TickLength', [0.03 0.01])
end
