clc
clear

a1 = 0.1:0.1:1;

for v=1:length(a1)
    a1(v)
    mi=0;
    all_periods=0;
    for g = 1:1000
        g
        file_path = 'as0 = %d_%0.5f_series.txt';
        sample = sprintf(file_path, g,a1(v));
        if ~exist(sample, 'file')
            warning('%s 不存在，跳过', sample);
            continue;
        end
        a = load(sample);
        t = a(:,1); x = a(:,2);
        x_smooth = smooth(x, 51, 'sgolay');
        [pks, locs] = findpeaks(x_smooth, t);
        rows=10

        peak_time(g,:)=locs(1:1:rows,1);
        mi=mi+locs(1:rows,1);

        tn = peak_time(g, :);
        periods = diff(tn);
        all_periods=mean(periods)+all_periods;
    end
    T = all_periods/(g)
    mi=(mi/g)';
    sigma2=zeros(rows,1)
    for k=1:rows
        k
        sigma2(k,1) =sum((peak_time(:, k) - mi(1,k)).^2);
        k=k+1;
    end
    sigma2=sigma2/(g-1);
    h=figure(v)
    plot(mi,sigma2,"kx",'LineWidth', 1, 'MarkerSize', 10)
    hold on
    p = polyfit(mi, sigma2, 1);
    slope = p(1);
    intercept = p(2);
    fprintf('拟合直线斜率为：%f\n', slope);
    hold on
    x_fit = linspace(min(mi), max(mi), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    hold off
    xlim([-10, 150])
    ylim padded
    ylabel('Peak variance, \sigma^2')
    xlabel('Time')

    title(['\alpha_{S0} = ', num2str(a1(v))])

    set(gca,'xtick',0:30:150)
    set(gca,'LineWidth',1.2,'Fontsize',24)
    set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
    set(gca,'XTickLabelRotation',0);
    set(gca,'YTickLabelRotation',0);

    filename = ['Peak diffusion constant_',  'a1=', num2str(a1(v), '%.5f'), '.pdf'];
    print(h, '-r600', '-dpdf', filename);

    D=slope;
    D_phi=(2*(pi^2)*D)/(T^2);
    result1(v,:)=[a1(v),D,D_phi];

end

h1=figure(length(a1)+1)

plot(result1(:,1),result1(:,2),"k-o",'LineWidth', 1, 'MarkerSize', 10)

plot(result1(:,1),result1(:,3),"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)

xlim([0,1.2])
ylabel('Phase diffusion')
xlabel('\alpha_{S0}')
set(gca,'xtick',0:0.2:1.2)
% set(gca,'ytick',0:0.5:1)
set(gca,'LineWidth',1.2,'Fontsize',24)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);
set(gca,'YTickLabelRotation',0);
ax = gca();
ax.YRuler.Exponent = -2;
print(h1, '-r600', '-dpdf', 'Phase diffusion constant');


