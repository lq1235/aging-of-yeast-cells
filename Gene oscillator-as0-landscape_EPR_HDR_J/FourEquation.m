clc
clear
[t,x]=ode45(@fun,[0,300],[10.9630341446514	112.165963264298	252.822731958915	135.751824484358]);%ode23s

figure(1)

% plot(t,x(:,5))
% hold on
% plot(t,x(:,1))
% hold on 
% plot(t,x(:,2))
% hold on
plot(t,x(:,3))
hold on
plot(t,x(:,4))
hold on

hold on
legend('S','H')
xlabel('t')
ylabel('x')



figure(2)
%相空间轨迹
plot(x(:,4),x(:,3))
xlim([0 1000]);
ylim([0 1200]);

xlabel('H')
ylabel('S')


function dx=fun(t,x)

  
%参数设置
as = 30.5;
ah = 183;
ah0 = 0.1;
as0 = 0.1;
beta = 3.7;%4.6
dm = 0.3;
dh = 3.8;
ds = 0.2;
Kh = 326;
Ks = 185;
n1 = 3;
n2 = 4.8;

mS=x(1);
mH=x(2);
S=x(3);
H=x(4);

 dx1dt = as0+as*H^n1/(Kh^n1+H^n1)-dm*mS;

  dx2dt = ah0+ah*Ks^n2/(Ks^n2+S^n2)-dm*mH;

    dx3dt = beta*mS-ds*S;

      dx4dt = beta*mH-dh*H;
    
  dx=[dx1dt;dx2dt;dx3dt;dx4dt];

end