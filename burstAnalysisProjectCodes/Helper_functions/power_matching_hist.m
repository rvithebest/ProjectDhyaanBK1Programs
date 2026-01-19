function [matched_sg_indices,matched_fg_indices]=power_matching_hist(power_gatherer_sg,power_gatherer_fg)
    edges=min([power_gatherer_sg, power_gatherer_fg]):0.5:max([power_gatherer_sg, power_gatherer_fg]);
    hist_sg=histcounts(power_gatherer_sg,edges);
    hist_fg=histcounts(power_gatherer_fg,edges);
    % Ensure equal number of electrodes in both classes by matching histograms
    min_hist = min(hist_sg, hist_fg);
    matched_sg_indices = [];
    matched_fg_indices = [];
    for k = 1:length(min_hist)
        if min_hist(k) > 0
            sg_indices_in_bin = find(power_gatherer_sg >= edges(k) & power_gatherer_sg < edges(k+1));
            fg_indices_in_bin = find(power_gatherer_fg >= edges(k) & power_gatherer_fg < edges(k+1));
    
            if length(sg_indices_in_bin) > min_hist(k)
                sg_indices_in_bin = randsample(sg_indices_in_bin, min_hist(k));
            end
            if length(fg_indices_in_bin) > min_hist(k)
                fg_indices_in_bin = randsample(fg_indices_in_bin, min_hist(k));
            end
    
            matched_sg_indices = [matched_sg_indices, sg_indices_in_bin];
            matched_fg_indices = [matched_fg_indices, fg_indices_in_bin];
        end
    end
end