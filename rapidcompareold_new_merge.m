%comp

for i = 1:321
    i
    
    diff = DMERGEnew(i,:) - DMERGE(i,:);
    different = diff >= 0.00000001  | diff <= -0.00000001;
    sum_diff = sum(different)
    
    miss_new = isnan(DMERGEnew(i,:));
    miss_old = isnan(DMERGE(i,:));
    
    diff_miss = sum(miss_new(1,:) - miss_old(1,:))
    
    
    if sum_diff ~= 0 || diff_miss ~= 0
        
        figure(1)
        plot(DMERGEnew(1,:), DMERGEnew(i,:), 'k.')
        hold on
        plot(DMERGE(1,:), DMERGE(i,:), 'ro')

        figure(2)
        plot(DMERGEnew(i,:), DMERGE(i,:), 'k.')
        good=~isnan(DMERGEnew(i,:)) & ~isnan(DMERGE(i,:));
        P=polyfit(DMERGEnew(i,good), DMERGE(i,good),1)

        figure(3)
        plot(DMERGE(1,:), DMERGEnew(i,:) -DMERGE(i,:), 'k.')

        close all
    end
    
end
