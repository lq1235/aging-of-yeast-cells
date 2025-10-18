clc;
clear;
close all;
%Xm=1.0;
% index1=3;
% index2=4;
dif = [0.5; 0.5];
ms = 200;
xmin = [0;0];
xmax = [500;500];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);

[X,Y] = meshgrid(x,y);

a1 =0:0.5:25;

num = length(a1);
file_path = 'pp_D1=%0.5f.txt';
for j = 1:num
h=figure(j);
sample=sprintf(file_path,a1(j));
% h=figure;
% S=2
% aij=2
% sample = 'pp_a12=1.50000.txt';
px = load(sample);
p = reshape(px(:,3),ms,ms);
sum(sum(p))
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;      
%第一种处理方法
PP = eq(Pi,0)+Pi;     
P_eps=min(min(PP));  
P = P_eps*eq(Pi,0)+Pi;
% %第二种方法
% eps=1.0e-0;
% P=Pi+eps;
P1=P/max(max(P));
U = -log(P1);
% U_max =25;%U_max =8.6;
% U=U.*(U<U_max)+U_max.*(U>U_max);
% U(U==U_max)=NaN;
surf(X,Y,U)
% mesh(X,Y,U);
% pcolor(X,Y,U)
shading interp;
colormap([jet(256)]);
hold on;
% caxis([-10,25]);%caxis([-2,8.6]);

height=-10

zlim([height, 20])
% hold on

% pcolor(X, Y, U);
% shading interp;
% colormap(jet(256));


% 使用 pcolor 绘制图像并设置 ZData 为 -2
h_pcolor = pcolor(X, Y, U);
shading interp;
colormap(jet(256));

% 设置 ZData 使得 pcolor 绘制在 z = -2 处
h_pcolor.ZData = height* ones(size(U)); % 将 pcolor 图层设置在 z = -2 处

% 添加 ColorBar
colorbarHandle = colorbar; % 获取 ColorBar 对象句柄
cb=colorbar;

setCMapRatio([1/22,2/22,4/22],[.1,.5,1])

% 调整 ColorBar 属性
cb.FontSize =20; % 设置字体大小
cb.LineWidth = 2; % 设置线宽

% caxis([min(min(U)),max(max(U))-6.5]);
set(gcf,'Position',[100 100 600 600]);
%这段代码会将当前图形窗口的位置设置在屏幕左下角 (100, 100) 的位置，并将窗口的大小设置为宽 600 像素、高 600 像素的正方形

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
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'YTickLabelRotation',0);%46是字体的旋转角度

% 绘制归一化流线
mg = 1:15:200;
ng = mg;

% 计算流线、梯度力线
dx = x(2) - x(1);
dy = y(2) - y(1);
[GPx, GPy] = gradient(P1, dx, dy);
Jx = FPx .* P1 - dif(1) * GPx;
Jy = FPy .* P1 - dif(2) * GPy;
E = Jy.^2 + Jx.^2;
JJx = Jx ./ (sqrt(E) + eps);
JJy = Jy ./ (sqrt(E) + eps);
Fx= dif(1)*GPx./P1;
Fy=dif(2)*GPy./P1;
F=Fx.^2+Fy.^2;
FFx=Fx./(sqrt(F)+eps);
FFy=Fy./(sqrt(F)+eps);

% 使用 quiver3 绘制流线，将 z = -5
quiver3(X(mg,ng), Y(mg,ng), height*ones(size(X(mg,ng))), FFx(mg,ng), FFy(mg,ng), zeros(size(X(mg,ng))), 0.5, 'color', 'k', 'LineWidth', 1);
hold on;
quiver3(X(mg,ng), Y(mg,ng), height*ones(size(X(mg,ng))), JJx(mg,ng), JJy(mg,ng), zeros(size(X(mg,ng))), 0.5, 'color', 'w', 'LineWidth', 1);

% 生成文件名
filename = ['Landscape_S=',  'D=', num2str(a1(j), '%.5f'), '.pdf'];
% filename = ['Landscape_S=', num2str(S), '_a12=', num2str(a1(j), '%.3f'), '.pdf'];
% 打印图像到文件
print(h, '-r600', '-dpdf', filename);
end


function setCMapRatio(varargin)
% @author:slandarer
if nargin==2
    ax=gca;
    oriRatio=sort(varargin{1});
    breakPnt=sort(varargin{2});
elseif nargin==3
    ax=varargin{1};
    oriRatio=sort(varargin{2});
    breakPnt=sort(varargin{3});
end
% 原始数据处理
try
    CLimit=get(ax,'CLim');
catch
end
try
    CLimit=get(ax,'ColorLimits');
catch
end
breakPnt=[CLimit(1),breakPnt,CLimit(2)];
newRatio=diff(breakPnt);
oriCMap=colormap(ax);
CLen=size(oriCMap,1);
newRatio=newRatio./diff([0,oriRatio,1]);
newRatio=round(newRatio./max(newRatio).*400);
oriRatio=[oriRatio,1];
% 最开始部分颜色条构造
tempCMap=oriCMap(1:ceil(oriRatio(1).*CLen),:);
CInd2=kron((1:size(tempCMap,1)-1)',ones(newRatio(1),1));
newCMap=tempCMap(CInd2,:);
CInd3=oriRatio(1).*CLen-size(tempCMap,1)+1;
CInd3=round(CInd3.*newRatio(1));
newCMap=[newCMap;repmat(tempCMap(end,:),[CInd3,1])];
% 循环添加新的颜色
for i=2:length(oriRatio)
    CInd1=round(newRatio(i).*(ceil(oriRatio(i-1).*CLen)-oriRatio(i-1).*CLen));
    if abs(ceil(oriRatio(i).*CLen)-oriRatio(i).*CLen)>0
        CInd2=ceil(oriRatio(i-1).*CLen)+1:ceil(oriRatio(i).*CLen)-1;
    else
        CInd2=ceil(oriRatio(i-1).*CLen)+1:ceil(oriRatio(i).*CLen);
    end
    CInd2=kron(CInd2',ones(newRatio(i),1));
    CInd3=round(newRatio(i).*(oriRatio(i).*CLen-floor(oriRatio(i).*CLen)));
    if ceil(oriRatio(i).*CLen)==ceil(oriRatio(i-1).*CLen)
        CInd1=[];
        CInd3=round(newRatio(i).*(oriRatio(i).*CLen-oriRatio(i-1).*CLen));
    end
    newCMap=[newCMap;
        repmat(oriCMap(ceil(oriRatio(i-1).*CLen),:),[CInd1,1]);
        oriCMap(CInd2,:);
        repmat(oriCMap(ceil(oriRatio(i).*CLen),:),[CInd3,1])];
end
colormap(ax,newCMap);
end

