function histogram_all_burst(length_all_elec_med,length_all_elec_cont,bin_num)
    % bin_num=10 for burst lengths
    % bin_num=7 for time center
    % 2 samplt t-test for burst length
    [h_length_all,p_length_all]=kstest2(length_all_elec_med,length_all_elec_cont);
    [count,edges]=histcounts(length_all_elec_med,bin_num,'Normalization','probability');
    [count2,edges2]=histcounts(length_all_elec_cont,bin_num,'Normalization','probability');
    bin_centers_sg=(edges(1:end-1)+edges(2:end))/2;
    bin_centers_fg=(edges2(1:end-1)+edges2(2:end))/2;
    plot(bin_centers_sg, count, '-', 'LineWidth', 2, 'Color', 'm', 'Marker', '*');
    hold on;
    plot(bin_centers_fg, count2, '-', 'LineWidth', 2, 'Color',[0.5 0 1], 'Marker', '*');
    xlabel('Burst lengths (s)')
    ylabel('Fraction of bursts')
    legend('Meditators','Controls')
    disp(p_length_all)
    % font size - 25 for the axes
    % set(gca,'FontSize',25)
    % display in text box above the plot- significance- format **** (p-value)
    % if p_length_all<0.0001
    %     text(0.5,0.5,['****(',num2str(p_length_all),')'],'FontSize',10)
    % elseif p_length_all<0.001
    %     text(0.5,0.5,['***(',num2str(p_length_all),')'],'FontSize',10)
    % elseif p_length_all<0.01
    %     text(0.5,0.5,['**(',num2str(p_length_all),')'],'FontSize',10)
    % elseif p_length_all<0.05
    %     text(0.5,0.5,['*(',num2str(p_length_all),')'],'FontSize',10)
    % else
    %     text(0.5,0.5,['N.S.(',num2str(p_length_all),')'],'FontSize',10)
    % end
end