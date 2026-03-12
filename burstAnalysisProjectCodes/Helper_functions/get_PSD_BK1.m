function [PSD_gatherer_med,PSD_gatherer_cont,freqVals_all]=get_PSD_BK1(curr_protocol)
    parent_data_folder= 'C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\Meditation_Project';
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyzing occipital electrodes for G2 protocol
    % protocol_list= [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 
    gridType='EEG';
    capType = 'actiCap64_UOL';
    [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
    occipital_elec=electrodeGroupList{strcmp(groupNameList,'Occipital')};
    psdVals_ST=cell(1,length(idx_list));
    psdVals_BL=cell(1,length(idx_list));
    freqVals_all=cell(1,length(idx_list));
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
        Fs=round(1/(timeVals(2)-timeVals(1)));
        goodTimePosST = find(timeVals>=st_range(1),1) + (1:round(Fs*diff(st_range)));
        goodTimePosBL = find(timeVals>=bl_range(1),1) + (1:round(Fs*diff(bl_range)));
        % Set up multitaper
        params.tapers=[1 1];
        params.pad= -1;
        params.Fs=Fs;
        params.fpass= [0 150];
        params.trialave=1;
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
            [psdValsST_temp,freqVals]=mtspectrumc(LFP_data_temp(:,goodTimePosST)',params);
            [psdValsBL_temp,~]=mtspectrumc(LFP_data_temp(:,goodTimePosBL)',params);
            psdVals_ST{i}=[psdVals_ST{i}; psdValsST_temp'];
            psdVals_BL{i}=[psdVals_BL{i}; psdValsBL_temp'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqVals_all{i}=freqVals;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pairedSubjectNameList = getPairedSubjectsBK1;
    meditator_list = pairedSubjectNameList(:,1);
    control_list = pairedSubjectNameList(:,2);
    PSD_gatherer_med=[];
    PSD_gatherer_cont=[];
    for i=1:length(meditator_list)
        med_name=meditator_list{i};
        cont_name=control_list{i};
        med_idx=find(strcmp(goodSubjectList,med_name));
        cont_idx=find(strcmp(goodSubjectList,cont_name));
        if isempty(med_idx) || isempty(cont_idx)
            continue;
        end
        if isempty(psdVals_ST{med_idx}) || isempty(psdVals_ST{cont_idx})
            % To avoid computation of bad subjects
            continue;
        end
        mean_psd_st_temp_med=mean(psdVals_ST{med_idx},1,'omitnan');
        mean_psd_bl_temp_med=mean(psdVals_BL{med_idx},1,'omitnan');
        mean_psd_st_temp_cont=mean(psdVals_ST{cont_idx},1,'omitnan');
        mean_psd_bl_temp_cont=mean(psdVals_BL{cont_idx},1,'omitnan');
        norm_psd_temp_med=10*(log10(mean_psd_st_temp_med)-log10(mean_psd_bl_temp_med));
        norm_psd_temp_cont=10*(log10(mean_psd_st_temp_cont)-log10(mean_psd_bl_temp_cont));
        PSD_gatherer_med=[PSD_gatherer_med; norm_psd_temp_med];
        PSD_gatherer_cont=[PSD_gatherer_cont; norm_psd_temp_cont];
    end
end