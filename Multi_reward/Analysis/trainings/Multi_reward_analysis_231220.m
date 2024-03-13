% Analysis of the behavioral experiment on 12/20/23
% There were 4 sessions. The timings of reward switching were detected
% visually and incorporated into the anlaysis in ad-hoc way.

%% Experiment info
n_session = 4;

exp_info.n_session   =   n_session;
exp_info.n_feeder    =   4;
reward_list = {"pear","banana"};
reward_session = [1,2,0,0; 2,1,0,0; 1,2,0,0; 1,2,0,0]; % row = session, col = flavor. number is corresponding to the flavor in reward_flavor. 0 = no flavor is allocated.
exp_info.reward_list = reward_list;
exp_info.reward     =   reward_session;

save('exp_info.mat','exp_info')

%% Directory
data_dir = pwd;

cortex_dir  =   dir(fullfile(data_dir,'cortex','*c3d'));
cortex_file =   fullfile(cortex_dir.folder,cortex_dir.name);
reward_dir  =   dir(fullfile(data_dir,'reward','*.mat'));
cdp_dir     =   dir(fullfile(data_dir,'cdp','*cdp*.txt'));
cdp_file    =   fullfile(cdp_dir.folder,cdp_dir.name);
serial_dir  =   dir(fullfile(data_dir,'*Serial_numbers*'));
serial_file =   fullfile(serial_dir.folder,serial_dir.name);

res_dir_path     =   fullfile(data_dir,'analysis');
res_dir     =   dir(res_dir_path);
if isempty(res_dir)
    mkdir(res_dir_path)
end

%% Extract cdp
cdp_extracted   =   dir(fullfile(res_dir_path,'extracted*cdp*'));
beh_dir         =   dir(fullfile(res_dir_path,'Ext_Behavior*','Extracted_Behavior*'));
if isempty(cdp_extracted)
    cd(cdp_dir.folder)
    serial_number = load(serial_file);
    serial_number = serial_number';
    ExtractCdp_AF_TY_v2('cdp',res_dir_path,'Include',serial_number) % Extract cdp data
    cd(res_dir_path)
    Process_extracted_cdp_TY_v0(res_dir_path)
    
    cd(data_dir)
    cdp_extracted   =   dir(fullfile(res_dir_path,'extracted*cdp*'));
    beh_dir         =   dir(fullfile(res_dir_path,'Ext_Behavior*','Extracted_Behavior*'));
    load(fullfile(cdp_extracted.folder,cdp_extracted.name))
    load(fullfile(beh_dir.folder,beh_dir.name))
    
    
elseif length(cdp_extracted) == 1
    disp('Behavior data are found in the directory.')
    load(fullfile(cdp_extracted.folder,cdp_extracted.name))
    load(fullfile(beh_dir.folder,beh_dir.name))
    disp('Behavior data are loaded.')
else
    error('Multiple files were contained in the directory. Check the file name.')
end

%% Synchronization and reward signals
analog_dir = dir(fullfile(cortex_dir.folder,'analog_signal.mat'));
if isempty(analog_dir)
    [~,~,AnalogSignals,AnalogFrameRate,~,~,~] = readC3D_analog(cortex_file); % read C3D data. Column 1 = reward, Column 2 = TTL, Column 3 = CDP server
    save(fullfile(cortex_dir.folder,'analog_signal.mat'),'AnalogSignals','AnalogFrameRate')
elseif length(analog_dir)==1
    load(fullfile(analog_dir.folder,analog_dir.name))
else
    error('Multiple files were contained in the directory. Check the file name.')
end

%% Detect 3sec TTL signals based on their onsets
% sync_3sec       =   AnalogSignals(:,2);
% 
% sync_onset      =   find(diff(sync_3sec) > 4); % onset of TTL
% sync_end        =   find(diff(sync_3sec) < -4); % end of TTL
% 
% % remove artifacts (usually caused by switching Master-9 on/off
% if length(sync_end) > length(sync_onset)
%     disp('Weird signals are detected in analog channel 2')
%     weird_TTL   =   find(sync_end(1:length(sync_onset))-sync_onset<0); % identify when the order of TTL onsets and ends is reversed
%     sync_onset  =   sync_onset(1:weird_TTL-1);
%     sync_end    =   sync_end(1:weird_TTL-1);
% elseif length(sync_onset) > length(sync_end)
%     disp('Weird signals are detected in analog channel 2')
%     weird_TTL   =   find(sync_end - sync_onset(1:length(sync_end))<0); % identify when the order of TTL onsets and ends is reversed
%     sync_onset  =   sync_onset(1:weird_TTL-1);
%     sync_end    =   sync_end(1:weird_TTL-1);
% end
% 
% sync_duration   =   sync_end - sync_onset; % interval should be 50 ms (= 6 samples)
% artifacts       =   find(sync_duration ~= median(sync_duration)); % identify when the TTL interval is deviated
% if length(artifacts) == 0
%     disp('TTL interval was stable across session.')
% elseif length(artifacts) == 1
%     disp('One artifact was detected.')
%     disp('Only the signals before the artifact was extracted.')
%     sync_onset      =   sync_onset(1:artifacts-1);
%     sync_end        =   sync_end(1:artifacts-1);
%     sync_duration   =   sync_duration(1:artifacts-1);
% elseif  length(artifacts) >= 2
%     error('Multiple artifacts were detected.')
% end
% 
% sync_interval   =   diff(sync_onset); % duration of each TTL signal
% 
% % FIGURE == Check TTL consistency
% figure
% tiledlayout(2,2,"TileSpacing","tight");
% % initial 1 min of TTL
% nexttile([1,2])
% plot((1/AnalogFrameRate:1/AnalogFrameRate:60),sync_3sec(1:AnalogFrameRate*60),'k') 
% title('TTL (3sec)')
% xlabel('Time (sec)')
% % Duration of each TTL
% nexttile
% histogram(sync_duration)
% xlabel('Duration (sec)')
% % Interval of TTL signals
% nexttile
% histogram(sync_interval)
% xlabel('TTL interval (sec)')


%% Pesudo-fibonacci TTL signals
% Align time based on the pseudo-Fibonacci TTLs
sync_fibonacci  =   AnalogSignals(:,3);
cortex_time     =   (1:length(sync_fibonacci))'/AnalogFrameRate; % time points by the cortex
% Detect cortex times corresponding to peaks in the sync signal separated by at least 2s
% [~, TTL_cortex_times] = findpeaks(normalize(sync_fibonacci,'zscore'),cortex_time,'MinPeakHeight',1,'MinPeakDistance',2);
% Detect network times corresponding to crossing 10% reference level in the sync signal
[~,TTL_cortex_times,~] =   risetime(normalize(sync_fibonacci,'zscore'),cortex_time);
if length(TTL_cortex_times)==length(true_TTL)
    TTL_cortex_times    =   TTL_cortex_times(true_TTL); % Avoid artifacts by turning the Master-9 off
end
rec_duration = 10800;   %approx rec duration (s), 10800 covers 3 hours
TTL_time_diff = [21; 13; 8; 5; 4];  %TTL delays in s
TTL_abs_times = [0; cumsum(repmat(TTL_time_diff,round(rec_duration*2/sum(TTL_time_diff)),1))];
TTL_abs_times = TTL_abs_times(1:length(TTL_cortex_times));

cortex_t = interp1(TTL_cortex_times,TTL_abs_times,cortex_time,'linear','extrap');   % interpolation of time based on the TTL signals
AnalogSignals(:,4)  =   cortex_t;

% Check the quality of time alignment by TTLs

% figure
% tiledlayout(3,1,'TileSpacing','tight')
% ax(1) = nexttile; plot(AnalogSignals(:,4),AnalogSignals(:,3))
% % xlim([-0.1,0.1])
% % ylim([-0.1,3.5])
% title('cortex')
% xlabel('Time (sec)')
% 
% ax(2) =nexttile; plot(sync_data(:,8),sync_data(:,3))
% % xlim([-0.1,0.1])
% % ylim([-0.1,2.1])
% title('cdp')
% xlabel('Time (sec)')
% 
% % shift between cdp and cortex time
% % [~, TTL_cortex_times2] = findpeaks(normalize(AnalogSignals(:,3),'zscore'),AnalogSignals(:,4),'MinPeakHeight',3,'MinPeakDistance',2); % cortex
% [~,TTL_cortex_times2,~] =   risetime(normalize(AnalogSignals(:,3),'zscore'),AnalogSignals(:,4)); % detected TTL by Cortex
% if length(TTL_cortex_times2)==length(true_TTL)
%     TTL_cortex_times2    =   TTL_cortex_times2(true_TTL);
% end
% % [~, TTL_cdp_times2] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,8),'MinPeakHeight',1,'MinPeakDistance',2,'MinPeakWidth',0.048); % cortex
% [~,TTL_cdp_times2,~]    =   risetime(normalize(sync_data(:,3),'zscore'),sync_data(:,8));
% if length(TTL_cdp_times2)==length(true_TTL)
%     TTL_cdp_times2    =   TTL_cdp_times2(true_TTL);
% end
% time_shift  =   TTL_cdp_times2 - TTL_cortex_times2;
% ax(3) = nexttile;
% plot(TTL_cortex_times2,time_shift)
% title('Peak shift between cortex and cdp')
% xlabel('# of TTL')
% ylabel('Lag: cdp - cortex (sec)')
% linkaxes(ax,'x')


% shift between cdp and cortex time
% [~, TTL_cortex_times2] = findpeaks(normalize(AnalogSignals(:,3),'zscore'),AnalogSignals(:,4),'MinPeakHeight',3,'MinPeakDistance',2); % cortex
[~,TTL_cortex_times2,~] =   risetime(normalize(AnalogSignals(:,3),'zscore'),AnalogSignals(:,4)); % detected TTL by Cortex
if length(TTL_cortex_times2)==length(true_TTL)
    TTL_cortex_times2    =   TTL_cortex_times2(true_TTL);
end
% [~, TTL_cdp_times2] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,8),'MinPeakHeight',1,'MinPeakDistance',2,'MinPeakWidth',0.048); % cortex
[~,TTL_cdp_times2,~]    =   risetime(normalize(sync_data(:,3),'zscore'),sync_data(:,8));
if length(TTL_cdp_times2)==length(true_TTL)
    TTL_cdp_times2    =   TTL_cdp_times2(true_TTL);
end
time_shift  =   TTL_cdp_times2 - TTL_cortex_times2;

figure
histogram(time_shift,'FaceColor','k')
title('CDP - Cortex')
xlabel('Time shift (sec)')
ylabel('Frequency')

% shift Ciholas Time by 1.05s (t = 0 corresponds to the first '3s' Master-9 TTL)  
AnalogSignals(:,4)  =   AnalogSignals(:,4)+1.05;    

%% Trouble shooting with TTL detections
% 
% 
% % Shift of peak from onset: cdp vs cortex
% figure
% t = tiledlayout(2,1);
% title(t,'Shift between detected peaks and peak onsets')
% 
% ax(1) = nexttile;
% [~, TTL_cdp_times2] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,8),'MinPeakHeight',1,'MinPeakDistance',2,'MinPeakWidth',0.048); % cortex
% TTL_cdp_onset   =   sync_data(diff(sync_data(:,3))>1,8);
% TTL_cdp_onset   =   TTL_cdp_onset(true_TTL);
% plot(TTL_cdp_times2-TTL_cdp_onset)
% title('cdp')
% xlabel('# of TTL')
% ylabel('Lag: peak - onset')
% 
% ax(2) = nexttile;
% cortex_time     =   (1:length(sync_fibonacci))'/AnalogFrameRate; % time points by the cortex
% % Detect cortex times corresponding to peaks in the sync signal separated by at least 2s
% [~, TTL_cortex_times2] = findpeaks(normalize(sync_fibonacci,'zscore'),cortex_time,'MinPeakHeight',1,'MinPeakDistance',2);
% TTL_cortex_times2   =   TTL_cortex_times2(true_TTL);
% TTL_cortex_onset    =   cortex_time(diff(sync_fibonacci)>1);
% TTL_cortex_onset    =   TTL_cortex_onset(true_TTL);
% plot(TTL_cortex_times2-TTL_cortex_onset)
% title('cortex')
% xlabel('# of TTL')
% ylabel('Lag: peak - onset')
% linkaxes(ax,'x')
% 
% % Test different method of detecting TTLs
% [~,sync_onset,~]   =   risetime(AnalogSignals(:,3),AnalogFrameRate);
% [~,sync_end,~]   =   falltime(AnalogSignals(:,3),AnalogFrameRate);
% sync_onset    =   sync_onset(true_TTL);
% sync_duration   =   sync_end - sync_onset;

%% Reward signals
% reward_delivery =   AnalogSignals(:,1);
% reward_t        =   AnalogSignals(:,4);
% 
% figure;
% findpeaks(normalize(reward_delivery,'zscore'),reward_t,'MinPeakHeight',3,'MinPeakDistance',0.15)
% [reward_vol,reward_times,reward_width,~]  =   findpeaks(normalize(reward_delivery,'zscore'),reward_t,'MinPeakHeight',3);

% Concatenate all reward signals.
feeder_data = [];
first_reward_ts_local = zeros(n_session,1);
for dd = 1:length(reward_dir)
    load(fullfile(reward_dir(dd).folder,reward_dir(dd).name),'bhv_data')
    n_trials = size(bhv_data.trials,1);
    feeder_data_field      =   bhv_data.fields;
    feeder_data_field{9} = 'bat_ID';
    feeder_data_field{14} = 'session';
    first_reward_ts_local(dd)     =   bhv_data.trials(find(bhv_data.trials(:,12),1,'first'),11); % the local timestamp of the first reward signal
    bhv_data.trials  =   [bhv_data.trials,repmat(dd,n_trials,1)]; % Session idx
    feeder_data = [feeder_data;bhv_data.trials];
end

% Detect the reward signals recorded in the cortex system. Notice that
% non-reward landings were not detected. These non-rewarding landings
% should be identified visually. This ad-hoc analysis is not optimal, so I
% may need to modify the recording script so that I do not need to open a
% new recording window when I switch smoothie syringes.
minpeakdistance = 1; % Initial value for the minimum peak distance. Decrease this value until all the peaks are detected
all_peak_detected = false; % The flag if all the peaks are detected.
while ~all_peak_detected
    [~, reward_times,w_reward_signal,~] = findpeaks(normalize(AnalogSignals(:,1),'zscore'),AnalogSignals(:,4),'MinPeakHeight',3,'MinPeakDistance',minpeakdistance); % Detection by findpeaks
    [~, reward_times2,~] =   risetime(normalize(AnalogSignals(:,1),'zscore'),AnalogSignals(:,4)); % detection by risetime
    % Validate the detection by comparing the # of detected peaks by two
    % different methods
    if length(reward_times) ~= length(reward_times2)
        minpeakdistance = minpeakdistance - 0.1; % Decrease the threshold by 0.1
    else
        fprintf('Final MinPeakDistance was %.1f \n',minpeakdistance)
        fprintf('%d peaks were detected. \n',length(reward_times))
        all_peak_detected = true;
    end
end

% Artifact detection
% Artifact pulse rised when new sessions were started. The width of the Detect the artifacts
all_artifacts_detected = false;
th  =   3;
while ~all_artifacts_detected
    reward_artifacts = w_reward_signal < mean(w_reward_signal) - th*std(w_reward_signal);
    n_detected_artifacts = sum(reward_artifacts);
    if n_detected_artifacts > n_session-1
        th = th+1; % increase threshold
    elseif n_detected_artifacts == n_session-1
        all_artifacts_detected = true; % end loops
        fprintf('%d artifacts were detected. It matches with the number of sessions. \n',n_detected_artifacts)
    else
        error('Some artifacts might be overlooked. Check the initial threshold value.')
    end
end

% Check if all the reward signals are detected/recorded
n_estimate_reward = length(reward_times) - n_detected_artifacts; % number of reward signals recorded by the Cortex
n_actual_reward   = sum(feeder_data(:,12)==1); % number of landings detected by beam

if n_estimate_reward == n_actual_reward
    disp('All rewarding landings were recorded and detected.')
else
    error('Some rewarding landings might be ignored. Check the raw data.')
end

% Add information to feeder data table
for trial = 1:size(feeder_data,1)
    feeder  =   feeder_data(trial,10);
    session =   feeder_data(trial,14);
    feeder_data(trial,15) = reward_session(session,feeder);
end
feeder_data_field{15} = 'reward_type';


% % Global timestamps assignment
% first_reward_ts_global  =   reward_times([true;reward_artifacts(1:end-1)]); % the first timestamps of true reward signals in each session
% time_offset_session = first_reward_ts_global - first_reward_ts_local; % The session start time in global time
% reward_times_corrected = reward_times(~reward_artifacts); % Exclude artifacts from detected reward signals
% for trial = 1:size(feeder_data,1)
%     session = feeder_data(trial,14);
%     feeder_data(trial,16) = round(feeder_data(trial,11) + time_offset_session(session),2);
% end
% feeder_data_field{16} = 'global_time';
% 
% figure
% estimated_reward_time = feeder_data(feeder_data(:,12)==1,16);
% histogram(estimated_reward_time - reward_times_corrected)
% xlabel('lag (sec): Estimated - Correct reward time')
% ylabel('Count')

% Assign global timestamp to reward signal
feeder_activation_ts_global = [];
% feeder signals in each session
for session = 1:n_session
    feeder_activation_ts_local_session = feeder_data(feeder_data(:,14)==session,11); % local timestamp of feeder activation
    feeder_activation_success_ts_local_session = feeder_data(feeder_data(:,14)==session & feeder_data(:,12)==1,11); % local timestamp of rewarded feeder activation
    
    artifact_idx = find(reward_artifacts); % idx of artifacts detected in analog signals
    if session == 1
        reward_ts_global_session = reward_times(1:artifact_idx(1)-1); % global timestamp of reward delivery
    elseif session == n_session
        reward_ts_global_session = reward_times(artifact_idx(session-1)+1:end); % global timestamp of reward delivery
    else
        reward_ts_global_session = reward_times(artifact_idx(session-1)+1:artifact_idx(session)-1); % global timestamp of reward delivery
    end
    feeder_activation_ts_global_session = interp1(feeder_activation_success_ts_local_session,reward_ts_global_session,feeder_activation_ts_local_session,'linear','extrap'); % interpolation of timestamps for non-rewarding feeder activation
    feeder_activation_ts_global = [feeder_activation_ts_global;feeder_activation_ts_global_session];
end
feeder_data(:,16) = round(feeder_activation_ts_global,2);
feeder_data_field{16} = 'global_time';

save(fullfile(res_dir_path,'feeder.mat'),'feeder_data_field','feeder_data','reward_list')



%% Detect all landings
f   =   cell(n_tags,10); % Save the extracted behavior at the feeders
f_fields = cell(1,10); % field name of f
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position

% Timing of landings was detected by Angelo's script
blanding = 1 - bflying;

% Compute the index of timestamps t for the start and end of landings for each bat.
% Also extract the 3D positions for each landing and compute the median of these positions.
for bb = 1:n_tags

    landing_start    =   find(diff([blanding(1,bb);blanding(:,bb)])==1);
    landing_end      =   find(diff([blanding(:,bb);blanding(end,bb)])==-1);

    if blanding(1,bb)~=1 % if the bat was flying when the recording started
        if blanding(end,bb)~=1 % if the bat was flying when the recording finished
        else % if the bat was landing when the recording finished
            landing_start(end) = [];
        end

    elseif blanding(1,bb)==1 % if the bat was landing when the recording started
        if blanding(end,bb)~=1 % if the bat was flying when the recording finished
            landing_end(1) = []; 
        else % if the bat was landing when the recording finished
            landing_start(end) = [];
            landing_end(1) = [];
        end
    end

    if length(landing_start) ~= length(landing_end)
        error('# of takeoff and landing mismatch. Check the detection.')
    else
        fprintf('%dthe bat: %d landing was detected. \n',bb, length(landing_start))
    end

    f{bb,1} = t(landing_start);
    f{bb,2} = t(landing_end);
    f_fields{1} = 'land_start_t';
    f_fields{2} = 'land_end_t';
    
    n_land = length(landing_start);
    f{bb,3} = n_land;
    f_fields{3} = 'n_land';
     
    r_land = cell(1,n_land);
    r_land_median = cell(1,n_land);
    for ll = 1:n_land
        r_land{ll}   =   r(landing_start(ll):landing_end(ll),:,:);
        r_land_median{ll}   =   median(r_land{ll}(:,:,bb));
    end
    f{bb,4} = r_land;
    f{bb,5} = r_land_median;
    f_fields{4} = 'r_land';
    f_fields{5} = 'r_land_med';

    % n_land    =   length(landing_start);
    % 
    % feeder_distance = NaN(n_land,3);
    % for ll = 1:n_land
    %     r_land   =   r(landing_start(ll):landing_end(ll),:,:);
    %     r_land_median   =   median(r_land(:,:,bb));
    % 
    %     [closest_feeder, feeder_id] = min(vecnorm(r_fd-r_land_median,2,2));
    %     feeder_distance(ll,1) = closest_feeder;
    %     feeder_distance(ll,2) = feeder_id;
    %     if closest_feeder < 0.5
    %         feeder_distance(ll,3) = 1;
    %     else
    %         feeder_distance(ll,3) = 0;
    %     end
    % end
    % f{bb,3} = feeder_distance;

    % nexttile;histogram(at_feeder{bb,3}(1,:),linspace(0,6,30))

    
    
    bfeeder = feeder_data(feeder_data(:,strcmp(feeder_data_field,'bat_ID'))==bb,[9,10,12,15,16]);
    bfeeder_fields = feeder_data_field([9,10,12,15,16]);

    
    f_fields{6} = 'activated_feeder';
    f{bb,6} = feeder_data(feeder_data(:,strcmp(feeder_data_field,'bat_ID'))==bb,10);
    f_fields{7} = 'outcome';
    f{bb,7} = feeder_data(feeder_data(:,strcmp(feeder_data_field,'bat_ID'))==bb,12);
    f_fields{8} = 'reward_type';
    f{bb,8} = feeder_data(feeder_data(:,strcmp(feeder_data_field,'bat_ID'))==bb,15);
    f_fields{9} = 'b_feeder_t';
    f{bb,9} = feeder_data(feeder_data(:,strcmp(feeder_data_field,'bat_ID'))==bb,16);
    

    % % Exclude the feeder activation during flights
    % t_feed = f{bb,6};
    % include_flag = true(length(t_feed),1);
    % for ll = 1:length(t_feed)
    %     if find(f{bb,1} <= t_feed(ll),1,'last') ~= find(f{bb,2} >= t_feed(ll),1,'first')
    %         include_flag(ll) = false;
    %     end
    % end
    % f{bb,6} = f{bb,6}(include_flag);
    % f{bb,7} = f{bb,7}(include_flag);
    % f{bb,8} = f{bb,8}(include_flag);
    % f{bb,9} = f{bb,9}(include_flag);

    % Median of landing positions when the feeder is activated
    t_feed = f{bb,9};
    r_land_median = f{bb,5};
    r_median_feed = cell(length(t_feed),1);
    for ll = 1:length(t_feed)
        r_median_feed{ll} = r_land_median{find(f{bb,1} <= t_feed(ll),1,'last')};
    end
    f{bb,10} = r_median_feed;
    f_fields{10} = 'r_median_feed';
end

% % Check if the feeders are activated while the bats land
% figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% tiledlayout(n_tags,1,"TileSpacing",'tight')
% 
% for bb = 1:n_tags
%     ax(bb) = nexttile;
%     area(t,bflying(:,bb)*5,'FaceAlpha',0.5,'LineStyle','none'); 
%     hold on
%     % area(t,ismember(round(t,2),round(f{bb,strcmp(f_fields,'b_feeder_t')},2))*3,'FaceAlpha',0.7,'LineStyle','none','FaceColor','r')
%     plot(t,r(:,1,bb))
%     scatter(round(f{bb,strcmp(f_fields,'b_feeder_t')},2),5*ones(length(f{bb,strcmp(f_fields,'b_feeder_t')}),1),10,'k')
%     hold off
%     legend({'Fly','x','feeder'})
% end
% xlabel('Tiem (sec)')
% linkaxes(ax)
% 
% % Coordinate of position when the feeders were activated
% figure();   set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.5]);
% tiledlayout(1,3,"TileSpacing",'tight')
% 
% for bb = 1:n_tags
%     r_med_feed = f{bb,10};
% 
%     nexttile(1); hold on
%     for ll = 1:size(r_med_feed,1)
%         scatter(r_med_feed{ll}(1),r_med_feed{ll}(2)) % XY
%     end
%     hold off
%     xlim(r_lim(1,:))
%     ylim(r_lim(2,:))
%     xlabel('x')
%     ylabel('y')
% 
%     nexttile(2); hold on
%     for ll = 1:size(r_med_feed,1)
%         scatter(r_med_feed{ll}(1),r_med_feed{ll}(3)) % XZ
%     end
%     hold off
%     xlim(r_lim(1,:))
%     ylim(r_lim(3,:))
%     xlabel('x')
%     ylabel('z')
% 
%     nexttile(3); hold on
%     for ll = 1:size(r_med_feed,1)
%         scatter(r_med_feed{ll}(2),r_med_feed{ll}(3)) % YZ
%     end
%     hold off
%     xlim(r_lim(2,:))
%     ylim(r_lim(3,:))
%     xlabel('y')
%     ylabel('z')
% 
%     pause(3)
% end

% %% Correct feeder positions (s)
% r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position
% 
% %=== Correct feeder position(s) 
% if Group_name == 'D'
%     [~,fd_idx] = min(vecnorm(r_fd-centroid_all,2,2));
%     r_fd(1,:) = centroid_all(fd_idx,:);                                                  
% else
%     for i=1:4
%     [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid_all,2,2));
%     if feeder_distance<0.2,r_fd(i,:) =  centroid_all(fd_idx,:);end  % Do not correct if further than 20 cm
%     end
% end

%% Visualize all landing
for bb = 1:n_tags

    t_feeder = f{bb,9};
    n_land_feeder = length(t_feeder);
    range_draw = 10; % range to visualize around the feeder activation timing
    
    figure('Units','normalized','OuterPosition',[0,0.05,1,0.9])
    tiledlayout(5,5,"TileSpacing",'loose')
    figure_count = 1;
    for ll = 1:n_land_feeder
        idx_feeder = find(ismember(round(t,2),round(t_feeder(ll),2)));
            if idx_feeder > 100*range_draw && idx_feeder < T - 100*range_draw
                r_plot = r(idx_feeder-100*range_draw:idx_feeder+100*range_draw,:,bb);
                t_plot = (-1*range_draw:0.01:range_draw)' + t(idx_feeder);
                l_plot = blanding(idx_feeder-100*range_draw:idx_feeder+100*range_draw,bb);
            elseif idx_feeder <= 100*range_draw
                r_plot = r(1:idx_feeder+100*range_draw,:,bb);
                t_plot = ((idx_feeder-1)/100:0.01:range_draw)' + t(idx_feeder);
                l_plot = blanding(1:idx_feeder+100*range_draw,bb);
            else
                r_plot = r(idx_feeder-100*range_draw:end,:,bb);
                t_plot = (-1*range_draw:0.01:(T-idx_feeder)/100)' + t(idx_feeder);
                l_plot = blanding(idx_feeder-100*range_draw:end,bb);
            end
        n = length(t_plot);
        cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    
        % t_plot = t(idx_feeder-100*60:idx_feeder+100*30);
    
        nexttile([1,3])
        hold on
        plot(t_plot,r_plot,'LineWidth',1);
        plot([t(idx_feeder),t(idx_feeder)],[-3,3],'k')
        area(t_plot,l_plot,'FaceAlpha',0.2,'LineStyle','none','FaceColor',[0.3010 0.7450 0.9330])
        hold off
        xlim([t(idx_feeder)-range_draw,t(idx_feeder)+range_draw])
        ylim([-3,3])
        xlabel('Time (sec)')
        legend({'x','y','z','Beam blocked'})
    
        nexttile([1,1])
        hold on
        p = plot(r_plot(:,1),r_plot(:,2));
        drawnow
        set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
        scatter(r(idx_feeder,1,bb),r(idx_feeder,2,bb),20,'k','filled')
        hold off
        xlabel('x')
        ylabel('y')
        xlim([r_lim(1,1),r_lim(1,2)])
        ylim([r_lim(2,1),r_lim(2,2)])
    
        nexttile([1,1])
        hold on
        p = plot(r_plot(:,2),r_plot(:,3));
        drawnow
        set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
        scatter(r(idx_feeder,2,bb),r(idx_feeder,3,bb),20,'k','filled')
        hold off
        xlabel('y')
        ylabel('z')
        xlim([r_lim(2,1),r_lim(2,2)])
        ylim([r_lim(3,1),r_lim(3,2)])
    
        if figure_count < 5
            figure_count = figure_count + 1;
        else
            if ll ~= n_land_feeder
                figure('Units','normalized','OuterPosition',[0,0.05,1,0.9])
                tiledlayout(5,5,"TileSpacing",'loose')
                figure_count = 1;
            end
        end
    end
    
    mkdir(fullfile(res_dir_path,'land_feeder'))
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        saveas(figHandles(i),fullfile(res_dir_path,['land_feeder_bat' num2str(bb) '_' num2str(i) '.png']));
    end
    close all;

end