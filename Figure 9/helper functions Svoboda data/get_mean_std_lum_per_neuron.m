function [meanlum, stdlum] =  get_mean_std_lum_per_neuron(lummat, plotyn)

[Nneuron, Ntime] = size(lummat);

if plotyn
    f = figure;
end
meanlum = nan(1,Nneuron);
stdlum = nan(1,Nneuron);
for nn=1:Nneuron
    if plotyn
        subplot(2,2,1)
        hold all
        plot(lummat(nn,:))
        xlabel('time step')
        ylabel('\Delta F / F')
        title('\Delta F / F over time')
        subplot(2,2,2)
        hold all
        [N,e] = histcounts(lummat(nn,:), 'Normalization','pdf');
        e = e(2:end) - (e(2)-e(1))/2;
        plot(e, N);    
        xlabel('\Delta F / F')
        ylabel('# measurements')
        title('Distribution \Delta F / F for each neuron')
    end
    meanlum(nn) = nanmean(lummat(nn,:));
    stdlum(nn) = nanstd(lummat(nn,:));
end
if plotyn
    subplot(2,2,3)
    Y= nanmean(lummat);
    x = 1:length(Y);
    Y_err = nanstd(lummat);                 
    errmat = [Y'-Y_err' 2*Y_err'];          %matrix for errorbars
    plot(Y, 'b', 'Linewidth', 2)                         %main plot
    hold on
    g=area(x, errmat);                                  %error bars
    set(g(1),'FaceColor', 'none', 'EdgeColor', 'none')  %make lower area-plot invisible
    set(g(2),'FaceColor', 'b', 'EdgeColor', 'none')     %colour for main error fill-in
    alpha(0.3);   %make transparent (alpha 1 makes opaque,0 is completely see-through)
    xlabel('time step')
    ylabel('\Delta F / F')
    title('Mean / std \Delta F / F over time')

    subplot(2,4,7)
    histogram(meanlum)
    title('Distribution (neurons) mean lum value')
    xlabel('Mean \Delta F / F value')

    subplot(2,4,8)
    histogram(stdlum)
    title('Distribution (neurons) std lum value')
    xlabel('STD \Delta F / F value')
end    