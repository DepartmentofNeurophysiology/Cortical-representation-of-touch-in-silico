function [activitymat, differencemat, counts, countsdiff] = plot_single_trial(xloc, yloc, lummat, timebin0, edges)

    bins = edges+(edges(2)-edges(1))/2;
    bins = bins(1:end-1);
    pointsize = 20;
    
    Nneuron = length(xloc);
    
    xlimits = [10*floor(min(xloc)/10), 10*ceil(max(xloc)/10)];
    ylimits = [10*floor(min(yloc)/10), 10*ceil(max(yloc)/10)];
    for tn = -1:1
        subplot(2,3,tn+2)
        hold all
        try
            scatter(xloc, yloc, pointsize, lummat(:,timebin0+tn), 'filled');
        catch
            keyboard
        end
        xlim(xlimits)
        ylim(ylimits)
        c = colorbar;
        caxis([-1 1]) % best to have it symmetrical around 0
        colormap jet
        title(['Time touch ' num2str(tn)])
        c.Label.String ='\Delta F/ F';
        xlabel('cell location (\mu m)')
        ylabel('cell location (\mu m)')            
    end
    
    activitymat = lummat(:,timebin0+1);
    differencemat = lummat(:,timebin0+1) - lummat(:,timebin0-1);
    
    subplot(2,3,4)
    scatter(xloc, yloc, pointsize, differencemat, 'filled');       
    cd = colorbar;
    caxis([-1 1])
    colormap jet
    xlim(xlimits)
    ylim(ylimits)
    title('Difference in activity')
    
    subplot(2,3,5)
    h = histogram(lummat(:,timebin0+1), edges, 'Normalization', 'probability');
    counts=h.Values;
    title('Activity')
    xlabel('\Delta F / F')
    ylabel('fraction of cells')
    xlim([-1, 5 ])
    disp(['Mean activity = ' num2str(nanmean(lummat(:,timebin0+1))), '; Standard deviation = ' num2str(nanstd(lummat(:,timebin0+1)))])
    
    subplot(2,3,6)
    hdiff = histogram(lummat(:,timebin0+1)-lummat(:,timebin0-1), edges, 'Normalization', 'probability');
    countsdiff = hdiff.Values;
    xlim([-2, 4 ])
    title(['Difference Activity, ' num2str(Nneuron)  ' neurons'])
    xlabel('\Delta (\Delta F / F)')
    ylabel('fraction of cells')
return