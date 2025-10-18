clc
clear
dim=4;
D=0.5;
xmin=[0,0,0,0];
xmax=[500,500,500,500];
a1 = 0.1:0.1:1; 
b = length(a1);

for g=1:1000
g
for  j= 1:b
    j
    as01 = a1(j); 
% file_path = 'pp_D1=%0.10f.txt';
%      sample=sprintf(file_path,a1(j)); 
file_path = 'as0 = %d_%0.5f_series.txt';
     sample=sprintf(file_path,g,a1(j)); 
    as = 30.5;
    ah = 183;

    ah0 = 0.1; 
    as0 = as01;
    beta = 3.7; 
    dm = 0.3;
    dh = 3.8;
    ds = 0.2;
    Kh = 326;
    Ks = 185;
    n1 = 3;
    n2 = 4.8;

    mS=24.2676387075027;
    mH=161.244264460529;
    S=324.816463610110;
    H=206.962295436469;
        dt=0.1;
        t=0;
        i=1;
        t0=200;
        steps=t0/dt;

        N1=1;N2=1;N3=1;
        while t<=t0
            noise = sqrt(2*D)*randn(1,dim);

            mS= mS+dt*(as0 + as * H^n1 / (Kh^n1 + H^n1) - dm * mS)+sqrt(dt)*noise(1,1);

            mH= mH+dt*(ah0 + ah * Ks^n2 / (Ks^n2 + S^n2) - dm * mH)+sqrt(dt)*noise(1,2);

            S= S+dt*(beta * mS - ds * S)+sqrt(dt)*noise(1,3);

            H= H+dt*(beta * mH - dh * H)+sqrt(dt)*noise(1,3);

            if mS<xmin(1,1)
                mS=2*xmin(1,1)-mS;
            elseif mS>xmax(1,1)
                mS=2*xmax(1,1)-mS;
            else
                mS=mS;
            end

            if mH<xmin(1,2)
                mH=2*xmin(1,2)-mH;
            elseif mH>xmax(1,2)
                mH=2*xmax(1,2)-mH;
            else
                mH=mH;
            end

            if S<xmin(1,3)
                S=2*xmin(1,3)-S;
            elseif S>xmax(1,3)
                S=2*xmax(1,3)-S;
            else
                S=S;
            end

            if H<xmin(1,4)
                H=2*xmin(1,3)-H;
            elseif H>xmax(1,3)
                H=2*xmax(1,3)-H;
            else
                H=H;
            end
            t=t+dt;
            x(i,:)=[t,mS,mH,S,H];
            i=i+1;
        end
           dlmwrite(sample, [x(:,[1,4,5]),], 'delimiter', '\t');
end

end

lightBlue = [0.43, 0.33, 0.99];
lightred = [0.98, 0.34, 0.34];

figure(1)
plot(x(:,1),x(:,2), 'Color', lightred,'LineWidth', 1)
hold on

figure(2)
plot(x(:,1),x(:,3), 'Color', lightBlue,'LineWidth', 1)
hold on

figure(3)
plot(x(:,1),x(:,4),'Color', 'green','LineWidth', 1)
hold on

figure(4)
plot(x(:,1),x(:,5),'Color', 'green','LineWidth', 1)
hold on


figure(5)
plot(x(:,5),x(:,4),'Color', 'green','LineWidth', 1)
hold on
xlim([0 1000]);
ylim([0 1200]);
xlabel('H')
ylabel('S')
