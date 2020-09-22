function plot_average(xloc, yloc, activitymat, differencemat, edges)

pointsize = 20;

[Ntrial, Nneuron] = size(activitymat);

avactivity = nanmean(activitymat);
avdifference = nanmean(differencemat);

xlimits = [10*floor(min(xloc)/10), 10*ceil(max(xloc)/10)];
ylimits = [10*floor(min(yloc)/10), 10*ceil(max(yloc)/10)];

subplot(2,2,1)
scatter(xloc, yloc, pointsize, avactivity, 'filled');  
cd = colorbar;
caxis([-1 1]) % best to have it symmetrical around 0
colormap jet
xlim(xlimits)
ylim(ylimits)
title('Average activity at t+1')
xlabel('cell location (\mu m)')
ylabel('cell location (\mu m)')  
cd.Label.String ='Average \Delta F/ F';

subplot(2,2,2)
scatter(xloc, yloc, pointsize, avdifference, 'filled');  
cd = colorbar;
caxis([-1 1]) % best to have it symmetrical around 0
colormap jet
xlim(xlimits)
ylim(ylimits)
title('Average difference in activity')
xlabel('cell location (\mu m)')
ylabel('cell location (\mu m)')  
cd.Label.String ='Average difference \Delta F/ F';

subplot(2,2,3)
histogram(avactivity, edges, 'Normalization', 'probability');
title(['Average activity, ' num2str(Nneuron)  ' neurons; ' num2str(Ntrial) ' trials'])
xlabel('\Delta F / F')
ylabel('fraction of cells')
xlim([-1, 2 ])

subplot(2,2,4)
histogram(avdifference, edges, 'Normalization', 'probability');
title(['Average difference in activity'])
xlabel('\Delta (\Delta F / F)')
ylabel('fraction of cells')
xlim([-1, 2 ])