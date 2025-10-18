clc;
clear;
seed=1;rng(seed);
dim=2;
a1=0.1:0.2:4.7;

h=0.01;
dc=0.0005;
steps=7000;
iters=7000;
n_traj=5000;
n_varpara=length(a1);
k1=1;  i1=10;
series1=1; series2=2;
for k=1:n_varpara
    k
    i = 1;
    func=@(x)force_flower(x,a1(k));
    x0=load('fixpoint.txt');
    x0=x0(k,2:end);
    x0=x0';

    variance1=zeros(2,1);
    amplitude1=zeros(2,1);
    while i<=n_traj
        [t,x]=sode_Heue(func,h,x0,iters,dc);

        if any(x(1,:)>1.4)
            continue;
        end

        if any(x(1,:)<0.2)
            continue;
        end

        if any(x(2,:)>1.4)
            continue;
        end

        if any(x(2,:)<0.2)
            continue;
        end

        if k == k1 && i == i1
            x1=x';
            dlmwrite('series.txt', [x1(:,[series1,series2]),], 'delimiter', '\t');
        end
        variance = var(x,1, 2);
        variance1=variance+variance1;
        amplitude=max(x,[],2)-min(x,[],2);
        amplitude1=amplitude1+amplitude;

        traj_1=x(series1,:);
        [acf, lags]=autocorr(traj_1,'NumLags',steps-1);
        t0=lags(1:1:end)*h;
        index(k,i)=min(find(acf<0));
        tau_ik(k,i)=trapz(t(1:index(k,i)-1),acf(1:index(k,i)-1));
        tau=mean(tau_ik,2);

        traj_2=x(series2,:);
        [acf1, lags]=autocorr(traj_2,'NumLags',steps-1);
        t0=lags(1:1:end)*h;
        index1(k,i)=min(find(acf1<0));
        tau_ik1(k,i)=trapz(t(1:index1(k,i)-1),acf1(1:index1(k,i)-1));
        tau1=mean(tau_ik1,2);
        if k == k1 && i == i1
            dlmwrite('autocorrelation.txt', [t0',  acf'], 'delimiter', '\t');
            dlmwrite('autocorrelation1.txt', [t0',  acf1'], 'delimiter', '\t');
        end
        traj1=x([series1,series2],1:end);
        [xcf,lags1]=crosscorr(traj1(1,:),traj1(2,:),'NumLags',steps-1);
        CXY=xcf(steps:1:2*steps-1);
        CYX=xcf(steps:-1:1);
        deltaC=CXY-CYX;
        t1=lags1(steps:2*steps-1)*h;
        n=3;
        DeltaC(k,i)=sqrt(trapz( t1(1:n) ,deltaC(1:n).^2 )/(h*n));
        DeltaC1=mean(DeltaC,2);
        i=i+1;
        if k == k1 && i == i1
            dlmwrite('crosscorrelation.txt', [t1',  CXY', CYX'], 'delimiter', '\t');
            dlmwrite('crosscorrelation-difference.txt', [t1',  deltaC'], 'delimiter', '\t');
        end

    end
    data(k,:)=variance1 /n_traj;
    data1(k,:)=amplitude1/n_traj;
end
filename='CSD_DCC.mat';
save(filename);
dlmwrite('tau-tau1-DeltaC.txt', [a1', tau, tau1, DeltaC1], 'delimiter', '\t');

figure(1)
plot(a1,tau,"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlabel('Glucose%','FontSize',27);
ylabel('\tau','FontSize',27);
set(gca,'xtick',0:1:6)
set(gca,'ytick',0:0.5:6)
xlim([-0.1 6])
ylim([1.251911046992085,3.751911046992085])

figure(2)
plot(a1,tau1,"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlabel('Glucose%','FontSize',27);
ylabel('\tau','FontSize',27);
set(gca,'xtick',0:1:6)
set(gca,'ytick',0:0.5:6)
xlim([-0.1 6])
ylim([1.230308582216262,3.980308582216262])

figure(3)
plot(a1,DeltaC1,"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'xtick',0:1:5)
xlabel('Glucose%','FontSize',27);
ylabel('\DeltaCC','FontSize',27);
xlim([-0.1 6])
ylim([0.001777469024453,0.002503469024453])
set(gca,'xtick',0:1:6)

figure(4)
plot(x(1,:),x(2,:))
xlim([0 2.5])
ylim([0 2.5])


figure(5)
plot(a1,data(:,2),"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2)
hold on
xlabel('Glucose%','FontSize',27);
ylabel('variance','FontSize',27);

figure(6)
plot(a1,data1(:,2),"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 6.2)
hold on
xlabel('Glucose%','FontSize',27);
ylabel('amplitude','FontSize',27);

function f=force_flower(x,a1)
f=zeros(length(x),1);
S=x(1);
H=x(2);
D = a1;
gamma2 =0.6;
St = 1.2*((0.03^3 / (0.03^3 + D^3)) * 0.0515 - 0.02036 * D + 1.6207);
Ht = (0.04^3 / (0.04^3 + D^3)) * 0.126 - 0.06 * D + 1.999;
f(1) = 0.01 * (St - S) + (S^3 / (0.52^3 + S^3)) * (St - S) - (0.1 + H) * S;
f(2) = 0.01 * (Ht - H) + (H^3 / (0.62^3 + H^3)) * (Ht - H) - (0.3 + gamma2 * S) * H;
end

function [t,x] = sode_Heue(func,h,x0,steps,dc)
dim = length(x0);
t = 0:h:(steps-1)*h;
h_sqrt=sqrt(h);
gx = sqrt(2*dc);
x(:,1) = x0;
Xm=3.0;
%  Xm1=ones(dim,1)*Xm;
for j = 1:steps-1
    noise = h_sqrt*randn(dim,1);
    gxn = gx*noise;
    aux = h*feval(func,x(:,j))+gxn;
    fxh =feval(func,x(:,j)+aux);
    xn = x(:,j)+ 0.5*(aux+h*fxh+gxn);

    if any(xn<0)
        x(:,j+1) = abs(xn);
    elseif any(xn > Xm)
        indices = find(xn > Xm);
        xn(indices) = 2 * Xm - xn(indices);
        x(:,j+1) = xn;
    else
        x(:,j+1) = xn;
    end
end
end