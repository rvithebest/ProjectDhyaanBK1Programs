clc;clear;close all;
[median_length_med_G1, median_length_cont_G1, mean_amplitude_med_G1, mean_amplitude_cont_G1] = getBurst_amp_duration('G1');
[median_length_med_G2, median_length_cont_G2, mean_amplitude_med_G2, mean_amplitude_cont_G2] = getBurst_amp_duration('G2');
[median_length_med_M2, median_length_cont_M2, mean_amplitude_med_M2, mean_amplitude_cont_M2] = getBurst_amp_duration('M2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PSD_gatherer_med_G1,PSD_gatherer_cont_G1,freqVals_all_G1]=get_PSD_BK1('G1');
[PSD_gatherer_med_G2,PSD_gatherer_cont_G2,freqVals_all_G2]=get_PSD_BK1('G2');
[PSD_gatherer_med_M2,PSD_gatherer_cont_M2,freqVals_all_M2]=get_PSD_BK1('M2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Figure
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(3,3,[0.12 0.12 0.85 0.85],0.07,0.1,0);
subplot(plotHandles(1,1));
plot(freqVals_all_G1{1},mean(PSD_gatherer_med_G1,1,'omitnan'),'-','LineWidth',2,'Color','m');
hold on;
plot(freqVals_all_G1{1},mean(PSD_gatherer_cont_G1,1,'omitnan'),'-','LineWidth',2,'Color','b');
xlabel('Frequency (Hz)');
xlim([0 80]);
ylabel('Mean Power (dB)');
% Draw dotted vertical lines at 24 and 34 Hz-slow gamma range
xline(24,'--k','LineWidth',1.5);
xline(34,'--k','LineWidth',1.5);
legend('Meditators','Controls');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,2));
violin_swarm_plot_paired(mean_amplitude_med_G1,mean_amplitude_cont_G1,0,1);
ylabel("Mean Amplitude");
subplot(plotHandles(1,3));
scatter_plot_power_burst(mean_amplitude_med_G1,mean_amplitude_cont_G1,median_length_med_G1 ...
        ,median_length_cont_G1,[],[]);
xlabel("Mean Amplitude");
ylabel("Median Burst length (s)")
[r_med_G1,p_med_G1]=corr(mean_amplitude_med_G1',median_length_med_G1',"Type","Spearman");
[r_con_G1,p_con_G1]=corr(mean_amplitude_cont_G1',median_length_cont_G1',"Type","Spearman");
legend({['Med: ','\rho = ',num2str(round(r_med_G1,2)),', p = ',num2str(round(p_med_G1,2))],...
        ['Cont: ','\rho = ',num2str(round(r_con_G1,2)),', p = ',num2str(round(p_con_G1,2))]},'Location','best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(2,1));
plot(freqVals_all_G2{1},mean(PSD_gatherer_med_G2,1,'omitnan'),'-','LineWidth',2,'Color','m');
hold on;
plot(freqVals_all_G2{1},mean(PSD_gatherer_cont_G2,1,'omitnan'),'-','LineWidth',2,'Color','b');
xlabel('Frequency (Hz)');
xlim([0 80]);
ylabel('Mean Power (dB)');
% Draw dotted vertical lines at 24 and 34 Hz-slow gamma range
xline(24,'--k','LineWidth',1.5);
xline(34,'--k','LineWidth',1.5);    
legend('Meditators','Controls');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(2,2));
violin_swarm_plot_paired(mean_amplitude_med_G2,mean_amplitude_cont_G2,0,1);
ylabel("Mean Amplitude");
subplot(plotHandles(2,3));
scatter_plot_power_burst(mean_amplitude_med_G2,mean_amplitude_cont_G2,median_length_med_G2 ...
        ,median_length_cont_G2,[],[]);
xlabel("Mean Amplitude");
ylabel("Median Burst length (s)")
[r_med_G2,p_med_G2]=corr(mean_amplitude_med_G2',median_length_med_G2',"Type","Spearman");
[r_con_G2,p_con_G2]=corr(mean_amplitude_cont_G2',median_length_cont_G2',"Type","Spearman");
legend({['Med: ','\rho = ',num2str(round(r_med_G2,2)),', p = ',num2str(round(p_med_G2,2))],...
        ['Cont: ','\rho = ',num2str(round(r_con_G2,2)),', p = ',num2str(round(p_con_G2,2))]},'Location','best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(3,1));
plot(freqVals_all_M2{1},mean(PSD_gatherer_med_M2,1,'omitnan'),'-','LineWidth',2,'Color','m');
hold on;
plot(freqVals_all_M2{1},mean(PSD_gatherer_cont_M2,1,'omitnan'),'-','LineWidth',2,'Color','b');
xlabel('Frequency (Hz)');
xlim([0 80]);
ylabel('Mean Power (dB)');
% Draw dotted vertical lines at 20 and 34 Hz-slow gamma range
xline(20,'--k','LineWidth',1.5);
xline(34,'--k','LineWidth',1.5);
legend('Meditators','Controls');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(3,2));
violin_swarm_plot_paired(mean_amplitude_med_M2,mean_amplitude_cont_M2,0,1);
ylabel("Mean Amplitude");
subplot(plotHandles(3,3));
scatter_plot_power_burst(mean_amplitude_med_M2,mean_amplitude_cont_M2,median_length_med_M2 ...
        ,median_length_cont_M2,[],[]);
xlabel("Mean Amplitude");
ylabel("Median Burst length (s)")
[r_med_M2,p_med_M2]=corr(mean_amplitude_med_M2',median_length_med_M2',"Type","Spearman");
[r_con_M2,p_con_M2]=corr(mean_amplitude_cont_M2',median_length_cont_M2',"Type","Spearman");
legend({['Med: ','\rho = ',num2str(round(r_med_M2,2)),', p = ',num2str(round(p_med_M2,2))],...
        ['Cont: ','\rho = ',num2str(round(r_con_M2,2)),', p = ',num2str(round(p_con_M2,2))]},'Location','best');
set_axis_ticks_fontsize(plotHandles,18,15,1);
set_axis_ticks_fontsize(plotHandles,18,15,2);
set_axis_ticks_fontsize(plotHandles,18,15,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox',...
[0.05 0.72 0.1 0.05],...
'String',{'Pre Med (G1)'},...
'FontWeight','bold',...
'Rotation',90,...
'FontSize',16,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox',...
[0.05 0.42 0.12 0.05],...
'String',{'Post Med (G2)'},...
'FontWeight','bold',...
'Rotation',90,...
'FontSize',16,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annotation('textbox',...
[0.05 0.1 0.12 0.05],...
'String',{'During Med (M2)'},...
'FontWeight','bold',...
'Rotation',90,...
'FontSize',16,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);