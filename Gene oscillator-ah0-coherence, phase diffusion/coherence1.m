clc
clear

a1 = 0.1:0.1:1;
b = length(a1);
for  j= 1:b
    j
    ah01 = a1(j); 
coherence_total=0;
    for g = 1:1000 
        g
        file_path = 'as0 = %d_%0.5f_series.txt';
        sample = sprintf(file_path, g,a1(j));
a=load(sample);
a(:, 2);
a(:, 3);
hold on


set(gca, 'color', 'white');
x_center = mean(a(1:end,2));

y_center = mean(a(1:end,3));

c0 = x_center + y_center*1i; 



r1 = a(1:end,2:3);



rr = r1(:,1)+r1(:,2)*1i-c0;
ang = angle(rr);
phi = diff(ang);
phi = phi - 2*pi*(phi>pi);
phi = phi + 2*pi*(phi<-pi);
B = sum(abs(phi));
phi0 = phi.*(phi>0);
A = sum(phi0);
coherence =abs(2*A/B - 1);
coherence_total=coherence_total+coherence
end
coherence11(j,:)=[a1(j),coherence_total/g];
end

h1=figure(1)
plot(coherence11(:,1),coherence11(:,2),"b-o",'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
xlabel('\fontsize{24} \alpha_{S0}');
ylabel('\fontsize{24} Coherence');

% ax = gca();
% ax.YRuler.Exponent = -1;
% ax.XRuler.Exponent = -1;
% set(gca,'xtick',0.125:0.01:0.165)
xlim([0,1.2])
ylim([0.947510695223913,1.002787202377028])
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
% axis([0 0.00395,0 0.00025])
set(gca,'XTickLabelRotation',0);
set(gca,'xtick',0:0.2:1.2)
set(gca,'ytick',0:0.01:1.2)
print(h1, '-r600', '-dpdf', 'Coherence');
