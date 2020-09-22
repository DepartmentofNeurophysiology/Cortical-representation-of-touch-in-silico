function figure_9 = combine_figs(InputDataStruct, volume, layername)

trialid_l23 = input(['Give the trial id of the figure to open for ' layername ' (Volume ' num2str(volume) ')']);
f_exp_single_trial = openfig([InputDataStruct.savefolder 'Volume ' num2str(volume) ' trial id ' num2str(trialid_l23) '.fig']);
axes_e = findall(f_exp_single_trial,'type','axes');
t_m1_e = axes_e(6);
t_p1_e = axes_e(4);
diff_e = axes_e(3);
act_e_hist = axes_e(2);
diff_e_hist = axes_e(1);

f_sim_single_trial = openfig([InputDataStruct.savefolder 'single_trial_' layername '_vol' num2str(volume) '_sims.fig']);
axes_s = findall(f_sim_single_trial,'type','axes');
t_m1_s = axes_s(6);
t_p1_s = axes_s(4);
diff_s = axes_s(3);
act_s_hist = axes_s(2);
act_s_hist = get(act_s_hist, 'Children');
set(act_s_hist, 'FaceColor', 'm')
diff_s_hist = axes_s(1);
diff_s_hist = get(diff_s_hist, 'Children');
set(diff_s_hist, 'FaceColor', 'm')

f_temp = figure;
for i=1:8 
    hAxRef(i)=subplot(2,4,i);
    pos{i} = {hAxRef(i).Position};
end
close(f_temp)

figure_9 = figure('Name', ['Compare data sims ' layername ' volume ' num2str(volume)]);
t_m1_e_new = copyobj(t_m1_e,figure_9);
set(t_m1_e_new,{'position'},pos{1}')
t_p1_e_new = copyobj(t_p1_e,figure_9);
set(t_p1_e_new,{'position'},pos{2}'')
diff_e_new = copyobj(diff_e,figure_9);
set(diff_e_new,{'position'},pos{3}')
act_e_hist_new = copyobj(act_e_hist,figure_9);
set(act_e_hist_new,{'position'},pos{4}'')
ha_e = get(act_e_hist_new, 'Children');
set(ha_e, 'FaceAlpha',1)
diff_e_hist_new = copyobj(diff_e_hist,figure_9);
set(diff_e_hist_new,{'position'},pos{8}')
hd_e = get(diff_e_hist_new, 'Children');
set(hd_e, 'FaceAlpha',1)


t_m1_s_new = copyobj(t_m1_s,figure_9);
set(t_m1_s_new,{'position'},pos{5}')
t_p1_s_new = copyobj(t_p1_s,figure_9);
set(t_p1_s_new,{'position'},pos{6}')
diff_s_new = copyobj(diff_s,figure_9);
set(diff_s_new,{'position'},pos{7}')
ha_s = copyobj(act_s_hist,act_e_hist_new);
hd_s = copyobj(diff_s_hist,diff_e_hist_new);

for i = 1:8
    subplot(2,4,i)
    box on
    set(gca, 'FontSize',16)
    if i==4 || i==8
        grid on
        legend('experiment', 'simulation')
        xlim([-1,2])
        ylim([0, 0.3])
    else
        colormap jet
    end    
end

close(f_exp_single_trial) 
close(f_sim_single_trial)