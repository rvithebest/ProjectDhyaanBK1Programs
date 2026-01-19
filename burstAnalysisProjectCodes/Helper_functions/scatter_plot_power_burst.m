function  scatter_plot_power_burst(power_gatherer_med,power_gatherer_cont,median_gatherer_med ...
        ,median_gatherer_cont,matched_sg_indices,matched_fg_indices)
    matched_sg_power = power_gatherer_med(matched_sg_indices);
    matched_fg_power = power_gatherer_cont(matched_fg_indices);
    matched_sg_lengths = median_gatherer_med(matched_sg_indices);
    matched_fg_lengths = median_gatherer_cont(matched_fg_indices);
    scatter(power_gatherer_med, median_gatherer_med, 50, 'm', 'o', 'LineWidth', 1.5);
    hold on;
    scatter(power_gatherer_cont, median_gatherer_cont, 50, 'o', 'LineWidth', 1.5,'MarkerEdgeColor',[0.5 0 1]);
    xlabel('Power difference (dB)')
    ylabel('Burst length (s)')
    % legend('Meditators','Controls');
    scatter(matched_sg_power, matched_sg_lengths, 50, 'mo','filled');
    scatter(matched_fg_power, matched_fg_lengths, 50, 'o', 'filled', 'MarkerFaceColor', [0.5 0 1]);
    legend('Meditators','Controls','Meditators (matched)','Controls (matched)')
end