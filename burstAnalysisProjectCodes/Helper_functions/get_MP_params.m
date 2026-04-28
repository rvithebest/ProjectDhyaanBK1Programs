function [st_range,bl_range,displayFlag,num_iterations,dict_size,adapt_dict_param,...
sg_freq,thresholdFraction,min_burst_num,u_bound]=get_MP_params(curr_protocol)
    % Intializes the params for running MP algorithm
    st_range= [0.25 1.25]; % stimulus time 
    % st_range=[0.25 1.2]; % - aliter_2
    bl_range=[-1 0]; % baseline time
    % bl_range=[-0.95 0]; % - aliter_2
    displayFlag=0;
    num_iterations=80;
    dict_size=2500000;
    adapt_dict_param=0.9;
    u_bound=1; % 1
    % slow gamma frequency range based on avg PSD plots
    if strcmp(curr_protocol,'G1') || strcmp(curr_protocol,'G2')
        % for G1 and G2 protocols
        sg_freq= [24 34];  %[24 34]
        thresholdFraction=1.7; %1.7 
        min_burst_num=100; % Total Num trials=120
    elseif strcmp(curr_protocol,'M2')
        % for M2 protocol
        sg_freq=[20 34]; 
        thresholdFraction=1.8; 
        min_burst_num=300; % Total Num trials=360
    else 
        error("Burst data only analyzed for G1, G2 and M2 protocols (gratings were displayed only while these protocols were run).");
    end 
end