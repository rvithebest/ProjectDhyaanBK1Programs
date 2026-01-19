clc;clear;close all;
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,2,[0.08 0.08 0.9 0.9],0.09,0.07,0);
% parent_data_folder='N:\Projects\ProjectDhyaan\BK1';  % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData
parent_data_folder= 'N:\Students\Vignesh\Meditation_Project';
analysis_data_folder= 'N:\Students\Vignesh\Meditation_Project\BK_1_EEG_MP_burst_analysis'; % Analysis results will be saved here
goodSubjectList=getGoodSubjectsBK1;
[allSubjectNames,expDateList]= getDemographicDetails('BK1');
badEyeCondition='ep'; % 'ep' : eye position
badTrialVersion= 'v8';
idx_list=1:length(goodSubjectList);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intialization
goodSubjectDate=cell(1,length(idx_list));
for i=idx_list
    subjectName=goodSubjectList{i};
    goodSubjectDate{i}=expDateList{strcmp(allSubjectNames,subjectName)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burst analysis parameters
st_range= [0.25 1.25]; % stimulus time
bl_range=[-1 0]; % baseline time
displayFlag=0;
num_iterations=80;
dict_size=2500000;
adapt_dict_param=0.9;
sg_freq=[24 34]; % slow gamma frequency range
thresholdFraction=1.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzing occipital electrodes for G2 protocol
% protocol_list= [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 
curr_protocol='G2';
gridType='EEG';
capType = 'actiCap64_UOL';
[electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
occipital_elec=electrodeGroupList{strcmp(groupNameList,'Occipital')};
length_accumulator=cell(1,length(idx_list));
onset_time_accumulator=cell(1,length(idx_list));
power_accumulator=cell(1,length(idx_list));
for i=1:length(idx_list)
    subjectName=goodSubjectList{i};
    disp(['Analyzing for subject: ' subjectName]);
    expDate=goodSubjectDate{i};
    protocol_name=curr_protocol;
    [bad_trials, bad_elec]= getBadTrialsAndElectrodes(subjectName,expDate,protocol_name,parent_data_folder,badEyeCondition,badTrialVersion);
    elec_list=1:64;
    num_elec=length(elec_list);
    segment_data_folder= fullfile(parent_data_folder,'data','segmentedData_v2',subjectName,gridType,expDate,protocol_name,'segmentedData_v2');
    timing_file=fullfile(segment_data_folder,'LFP','lfpInfo.mat');
    if ~exist(timing_file,'file')
        continue;
    end
    t=load(timing_file);
    timeVals=t.timeVals;
    e=load(fullfile(segment_data_folder,'LFP','elec1.mat'));
    good_trials=setdiff(1:size(e.analogData,1),bad_trials);
    num_trials=length(good_trials);
    if (num_trials<=30)
        continue;
    end
    good_elec_num=length(setdiff(occipital_elec,bad_elec));
    if good_elec_num<3
        continue;
    end
    for k=1:num_elec
        elec_label=elec_list(k);
        if ismember(elec_label,bad_elec) || ~ismember(elec_label,occipital_elec)
            continue;
        end
        e=load(fullfile(segment_data_folder,'LFP',['elec' num2str(elec_label) '.mat']));
        LFP_data_temp=e.analogData(good_trials,:);
        diffPower=getChangeInPower(LFP_data_temp,timeVals,st_range,bl_range,sg_freq);
        power_accumulator{i}=[power_accumulator{i}, (diffPower)];
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        save_folder=fullfile(analysis_data_folder,subjectName,expDate,protocol_name,['Elec' num2str(elec_label)]);
        load(fullfile(save_folder,'MP_results.mat'));
        [length_temp_sg,~,time_center_temp_sg,gaborInfo,header,~]=getBurstLengthMP(LFP_data_temp,timeVals,thresholdFactor,displayFlag,st_range,bl_range,sg_freq,num_iterations,adapt_dict_param,dict_size,gaborInfo,header,i);
        length_temp_all_trials=[];
        onset_temp_all_trials=[];
        for ii=1:length(length_temp_sg)
            if isempty(length_temp_sg{ii})
                continue;
            end
            % reject burst lengths longer than 1.25 second (Stimulus Period)
            reject_idx=find((length_temp_sg{ii}')>1.25);
            length_temp_sg{ii}(reject_idx)=[];
            time_center_temp_sg{ii}(reject_idx)=[];
            length_temp_all_trials=[length_temp_all_trials, length_temp_sg{ii}'];
            % Find the first burst (in terms of time) in each trial
            % Onset time (earliest)
            onset_time_temp_trial_sg=(time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5));
            onset_time_temp_trial_sg(onset_time_temp_trial_sg<0)=[];
            onset_idx=(find((onset_time_temp_trial_sg)==min(onset_time_temp_trial_sg)));
            onset_temp_all_trials=[onset_temp_all_trials,onset_time_temp_trial_sg(onset_idx)];
        end
        length_accumulator{i}=[length_accumulator{i}, (length_temp_all_trials)];
        onset_time_accumulator{i}=[onset_time_accumulator{i}, (onset_temp_all_trials)];
    end
end
pairedSubjectNameList = getPairedSubjectsBK1;
meditator_list = pairedSubjectNameList(:,1);
control_list = pairedSubjectNameList(:,2);
median_length_med=zeros(1,length(meditator_list));
median_length_cont=zeros(1,length(control_list));
mean_power_med=zeros(1,length(meditator_list));
mean_power_cont=zeros(1,length(control_list));
mean_onset_med=zeros(1,length(meditator_list));
mean_onset_cont=zeros(1,length(control_list));
length_all_med=[];
length_all_cont=[];
for i=1:length(meditator_list)
    med_name=meditator_list{i};
    cont_name=control_list{i};
    med_idx=find(strcmp(goodSubjectList,med_name));
    cont_idx=find(strcmp(goodSubjectList,cont_name));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(med_idx) || isempty(cont_idx)
        continue;
    end
    if isempty(length_accumulator{med_idx}) || isempty(length_accumulator{cont_idx})
        % To omit computation of bad subjects
        continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if less than 100 bursts- we do not conisder the pair
    if (length(length_accumulator{med_idx})<100) || (length(length_accumulator{cont_idx})<100)
         continue;
    end
    median_length_med(i)=median(length_accumulator{med_idx});
    median_length_cont(i)=median(length_accumulator{cont_idx});
    length_all_med=[length_all_med, length_accumulator{med_idx}];
    length_all_cont=[length_all_cont, length_accumulator{cont_idx}];
    mean_power_med(i)=10*log10(mean((power_accumulator{med_idx})));
    mean_power_cont(i)=10*log10(mean((power_accumulator{cont_idx})));
    mean_onset_med(i)=mean(onset_time_accumulator{med_idx});
    mean_onset_cont(i)=mean(onset_time_accumulator{cont_idx});
end
% Remove zeros
median_length_med=median_length_med(median_length_med~=0);
median_length_cont=median_length_cont(median_length_cont~=0);
mean_power_med=mean_power_med(mean_power_med~=0);
mean_power_cont=mean_power_cont(mean_power_cont~=0);
mean_onset_med=mean_onset_med(mean_onset_med~=0);
mean_onset_cont=mean_onset_cont(mean_onset_cont~=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,1));
histogram_all_burst(length_all_med,length_all_cont,10);
subplot(plotHandles(1,2));
violin_swarm_plot_paired(median_length_med,median_length_cont,0,0);
% violin_swarm_plot_paired(mean_power_med,mean_power_cont,0,1);
% ylabel("Mean Power (dB)")
% violin_swarm_plot_paired(mean_onset_med,mean_onset_cont,1,1);
% ylabel("Mean Onset Time (s)")
[matched_med_indices,matched_cont_indices]=power_matching_hist(mean_power_med,mean_power_cont);
subplot(plotHandles(2,1));
scatter_plot_power_burst(mean_power_med,mean_power_cont,median_length_med ...
        ,median_length_cont,matched_med_indices,matched_cont_indices);
% No significant correlation for both the groups
[r_med,p_med]=corr(mean_power_med',median_length_med',"Type","Spearman");
[r_con,p_con]=corr(mean_power_cont',median_length_cont','Type','Spearman');
subplot(plotHandles(2,2));
violin_swarm_plot(median_length_med(matched_med_indices),median_length_cont(matched_cont_indices));
set_axis_ticks_fontsize(plotHandles,18,15,1);
set_axis_ticks_fontsize(plotHandles,18,15,2);
figure;
violin_swarm_plot_paired(mean_power_med,mean_power_cont,0,1);