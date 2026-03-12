function [median_length_med, median_length_cont, mean_amplitude_med, mean_amplitude_cont] = getBurst_amp_duration(curr_protocol)
    parent_data_folder= 'C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\Meditation_Project';
    analysis_data_folder= 'C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\Meditation_Project\BK_1_EEG_MP_burst_analysis'; % Analysis results will be saved here
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
    gridType='EEG';
    capType = 'actiCap64_UOL';
    elec_list=1:64;
    num_elec=length(elec_list);
    [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
    occipital_elec=electrodeGroupList{strcmp(groupNameList,'Occipital')};
    length_accumulator=cell(1,length(idx_list));
    amplitude_accumulator=cell(1,length(idx_list));
    % Burst analysis parameters
    [st_range,bl_range,displayFlag,num_iterations,dict_size,adapt_dict_param,...
    sg_freq,thresholdFraction,min_burst_num]=get_MP_params(curr_protocol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyzing occipital electrodes
    for i=idx_list
        subjectName=goodSubjectList{i};
        disp(['Analyzing for subject: ' subjectName]);
        expDate=goodSubjectDate{i};
        protocol_name=curr_protocol;
        [bad_trials, bad_elec]= getBadTrialsAndElectrodes(subjectName,expDate,protocol_name,parent_data_folder,badEyeCondition,badTrialVersion);
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
        if (good_elec_num<3)
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
            thresholdFactor=sqrt(thresholdFraction*diffPower);
            save_folder=fullfile(analysis_data_folder,subjectName,expDate,protocol_name,['Elec' num2str(elec_label)]);
            load(fullfile(save_folder,'MP_results.mat'));
            [length_temp_sg,~,~,gaborInfo,header,~,amp_temp_sg]=getBurstLengthMP_2(LFP_data_temp,timeVals,thresholdFactor,displayFlag,st_range,bl_range,sg_freq,num_iterations,adapt_dict_param,dict_size,gaborInfo,header,i);
            length_temp_all_trials=[];
            amp_temp_all_trials=[];
            for ii=1:length(length_temp_sg)
                if isempty(length_temp_sg{ii})
                    continue;
                end
                % rejecting burst lengths longer than 1 seconds
                reject_idx=find((length_temp_sg{ii}')>1);
                length_temp_sg{ii}(reject_idx)=[];
                amp_temp_sg{ii}(reject_idx)=[];
                length_temp_all_trials=[length_temp_all_trials, length_temp_sg{ii}'];
                amp_temp_all_trials=[amp_temp_all_trials, amp_temp_sg{ii}'];
                % Find the first burst (in terms of time) in each trial
            end
            length_accumulator{i}=[length_accumulator{i}, (length_temp_all_trials)];
            amplitude_accumulator{i}=[amplitude_accumulator{i}, (amp_temp_all_trials)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pairedSubjectNameList = getPairedSubjectsBK1;
    meditator_list = pairedSubjectNameList(:,1);
    control_list = pairedSubjectNameList(:,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    median_length_med=zeros(1,length(meditator_list));
    median_length_cont=zeros(1,length(control_list));
    mean_amplitude_med=zeros(1,length(meditator_list)); 
    mean_amplitude_cont=zeros(1,length(control_list));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        if (length(length_accumulator{med_idx})<min_burst_num) || (length(length_accumulator{cont_idx})< min_burst_num)
             continue;
        end
        median_length_med(i)=median(length_accumulator{med_idx});
        median_length_cont(i)=median(length_accumulator{cont_idx});
        mean_amplitude_med(i)=mean(amplitude_accumulator{med_idx});
        mean_amplitude_cont(i)=mean(amplitude_accumulator{cont_idx});
    end
    % Remove zeros
    median_length_med=median_length_med(median_length_med~=0);
    median_length_cont=median_length_cont(median_length_cont~=0);
    mean_amplitude_med=mean_amplitude_med(mean_amplitude_med~=0);
    mean_amplitude_cont=mean_amplitude_cont(mean_amplitude_cont~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end