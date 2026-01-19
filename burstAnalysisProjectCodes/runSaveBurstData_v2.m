clc;clear;close all;
% Parallel pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
end
pool = parpool("Processes");
% parent_data_folder='N:\Projects\ProjectDhyaan\BK1';  % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData
parent_data_folder= 'D:\Vignesh\Meditation_Project';
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
gaborInfo_accumulator=cell(length(idx_list));
header_accumulator=cell(length(idx_list));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burst analysis parameters
st_range= [0.25 1.25]; % stimulus time
bl_range=[-1 0]; % baseline time
displayFlag=0;
num_iterations=80;
dict_size=2500000;
adapt_dict_param=0.9;
sg_freq=[20 35]; % slow gamma frequency range
thresholdFraction=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzing occipital electrodes for G2 protocol
% protocol_list= [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 
curr_protocol='G2';
gridType='EEG';
capType = 'actiCap64_UOL';
[electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
occipital_elec=electrodeGroupList{strcmp(groupNameList,'Occipital')};
parfor i=idx_list
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
    if num_trials==0
        continue;
    end
    gaborInfo_elec_gatherer=cell(1,num_elec);
    header_elec_gatherer=cell(1,num_elec);
    for k=1:num_elec
        elec_label=elec_list(k);
        if ismember(elec_label,bad_elec) || ~ismember(elec_label,occipital_elec)
            continue;
        end
        save_folder=fullfile(analysis_data_folder,subjectName,expDate,protocol_name,['Elec' num2str(elec_label)]);
        save_file=fullfile(save_folder,'MP_results.mat');
        if exist(save_file,'file')
            continue;
        end
        e=load(fullfile(segment_data_folder,'LFP',['elec' num2str(elec_label) '.mat']));
        LFP_data_temp=e.analogData(good_trials,:);
        diffPower=getChangeInPower(LFP_data_temp,timeVals,st_range,bl_range,sg_freq);
        thresholdFactor=sqrt(thresholdFraction*diffPower);
        uniq_idx=[subjectName '_' expDate '_Elec' num2str(elec_label) '_' protocol_name];
        [~,~,~,gaborInfo,header,~]=getBurstLengthMP(LFP_data_temp,timeVals,thresholdFactor,displayFlag,st_range,bl_range,sg_freq,num_iterations,adapt_dict_param,dict_size,[],[],uniq_idx);
        gaborInfo_elec_gatherer{k}=gaborInfo;
        header_elec_gatherer{k}=header;
        save_data_parfor(save_folder,gaborInfo,header);
    end
    gaborInfo_accumulator{i}=gaborInfo_elec_gatherer;
    header_accumulator{i}=header_elec_gatherer;
end
delete(pool);
% Saving the results
function save_data_parfor(save_folder,gaborInfo,header)
    if ~exist(save_folder,'dir')
        mkdir(save_folder)
    end
    save_file=fullfile(save_folder,'MP_results.mat');
    if exist(save_file,'file')
        return;
    end
    save(save_file,'gaborInfo','header');
end