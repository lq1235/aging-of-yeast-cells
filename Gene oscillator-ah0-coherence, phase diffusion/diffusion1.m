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

          rows = size(locs, 1);  
        rows=10

        peak_time(g,:)=locs(1:1:rows,1);
        mi=mi+locs(1:rows,1);

        tn = peak_time(g, :);
        periods = diff(tn);  
        all_periods=mean(periods)+all_periods;
    end

h=figure(v)

for g1 = 1:size(peak_time, 1)

    tn = peak_time(g1, :);  
    scatter(tn, g1 * ones(size(tn)), 10, 'filled', 'k')  
    hold on
end

xlabel('Peak time')
ylabel('Trajectory index')
title(['\alpha_{S0} = ', num2str(a1(v))])
% xlim([min(peak_time(:)), max(peak_time(:))])
xlim([-10, 150])
ylim([0 1000])
ax = gca;
ax.Box = 'on'; 
ax.XAxisLocation = 'bottom'; 
ax.YAxisLocation = 'left'; 

set(gca,'xtick',0:30:150)
set(gca,'LineWidth',1.2,'Fontsize',24)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);
set(gca,'YTickLabelRotation',0);

filename = ['Raster plot of the peak times_',  'a1=', num2str(a1(v), '%.5f'), '.pdf'];
print(h, '-r600', '-dpdf', filename);

end





