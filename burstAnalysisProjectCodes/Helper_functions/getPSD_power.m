function [norm_power_PSD,peak_freq_PSD,psdVals_ST,psdVals_BL,freqVals_all]=getPSD_power(curr_protocol)
    % parent_data_folder='N:\Projects\ProjectDhyaan\BK1';  % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData
    % parent_data_folder= 'N:\Students\Vignesh\Meditation_Project';
    % analysis_data_folder= 'N:\Students\Vignesh\Meditation_Project\BK_1_EEG_MP_burst_analysis'; % Analysis results will be saved here
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Burst analysis parameters
    st_range= [0.25 1.25]; % stimulus time
    bl_range=[-1 0]; % baseline time
    sg_freq=[24 34]; % slow gamma frequency range
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyzing occipital electrodes for G2 protocol
    % protocol_list= [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 
    gridType='EEG';
    capType = 'actiCap64_UOL';
    [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
    occipital_elec=electrodeGroupList{strcmp(groupNameList,'Occipital')};
    psdVals_ST=cell(1,length(idx_list));
    psdVals_BL=cell(1,length(idx_list));
    peak_freqVals=cell(1,length(idx_list));
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
            sg_freq_pos= intersect(find(freqVals>=sg_freq(1)),find(freqVals<sg_freq(2)));
            norm_psd_temp=10*log10(psdValsST_temp(sg_freq_pos)./psdValsBL_temp(sg_freq_pos));
            norm_psd_temp=norm_psd_temp';
            [~,peak_freq_pos]=max(norm_psd_temp,[],2);
            slow_gamma_freqVals=freqVals(sg_freq_pos);
            peak_freqVals_temp=slow_gamma_freqVals(peak_freq_pos);
            if (peak_freqVals_temp==sg_freq(1)) || (peak_freqVals_temp==sg_freq(2))
                continue;
            end
            peak_freqVals{i}=[peak_freqVals{i}; peak_freqVals_temp];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    norm_power_PSD=cell(1,length(idx_list));
    peak_freq_PSD=zeros(1,length(idx_list));
    for i=1:length(idx_list)
        if isempty(psdVals_ST{i})
            continue;
        end
        mean_PSD_st=mean(psdVals_ST{i},1,'omitnan');
        mean_PSD_bl=mean(psdVals_BL{i},1,'omitnan');
        freq_pos= intersect(find(freqVals_all{i}>=sg_freq(1)),find(freqVals_all{i}<sg_freq(2)));
        log10_PSD_st=log10(sum(mean_PSD_st(:,freq_pos),2));
        log10_PSD_bl=log10(sum(mean_PSD_bl(:,freq_pos),2));
        norm_power_PSD{i}=10*(log10_PSD_st - log10_PSD_bl);
        % Peak frequency
        % mean_PSD_norm=mean_PSD_st./mean_PSD_bl;
        % [~,peak_freq_pos]=max(mean_PSD_norm(freq_pos));
        % sg_freq_vals=freqVals_all{i}(freq_pos);
        % peak_freq_PSD(i)=sg_freq_vals(peak_freq_pos);
        peak_freq_PSD(i)=mean(peak_freqVals{i});
    end
end