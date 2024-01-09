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

res_dir_path     =   fullfile(data_dir,'analysis');
res_dir     =   dir(res_dir_path);
if isempty(res_dir)
    mkdir(res_dir_path)
end

%% Extract cdp
cdp_extracted   =   dir(fullfile(res_dir_path,'extracted*'));
beh_dir         =   dir(fullfile(res_dir_path,'Ext_Behavior*','Extracted_Behavior*'));
if isempty(cdp_extracted)
    cd(cdp_dir.folder)
    ExtractCdp_AF_TY_v2('cdp',res_dir_path) % Extract cdp data
    Process_extracted_cdp_TY_v0(res_dir_path)
    cdp_extracted   =   dir(fullfile(cdp_dir.folder,'extracted*'));
    beh_dir         =   dir(fullfile(cdp_dir.folder,'Ext_Behavior*','Extracted_Behavior*'));
    load(fullfile(cdp_extracted.folder,cdp_extracted.name))
    load(fullfile(beh_dir.folder,beh_dir.name))
    cd(data_dir)
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
    feeder_data_field      =   [bhv_data.fields,'session'];
    first_reward_ts_local(dd)     =   bhv_data.trials(find(bhv_data.trials(:,12),1,'first'),11); % the local timestamp of the first reward signal
    bhv_data.trials  =   [bhv_data.trials,repmat(dd,n_trials,1)];
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


% Global timestamps assignment
first_reward_ts_global  =   reward_times([true;reward_artifacts(1:end-1)]); % the first timestamps of true reward signals in each session
time_shift_session = first_reward_ts_global - first_reward_ts_local; % The session start time in global time
reward_times_corrected = reward_times(~reward_artifacts); % Exclude artifacts from detected reward signals
for trial = 1:size(feeder_data,1)
    session = feeder_data(trial,14);
    feeder_data(trial,15) = round(feeder_data(trial,11) + time_shift_session(session),2);
end
feeder_data_field{15} = 'global_time';

figure
estimated_reward_time = feeder_data(feeder_data(:,12)==1,15);
histogram(estimated_reward_time - reward_times_corrected)
xlabel('lag (sec): Estimated - Correct reward time')
ylabel('Count')

% Add information to feeder data table
for trial = 1:size(feeder_data,1)
    feeder  =   feeder_data(trial,10);
    session =   feeder_data(trial,14);
    feeder_data(trial,16) = reward_session(session,feeder);
end
feeder_data_field{16} = 'reward_type';

save(fullfile(res_dir_path,'extracted_feeder.mat'),'feeder_data_field','feeder_data','reward_list')

%% Detect the flights to the feeders
at_feeder   =   cell(n_tags,10);
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position

blanding = 1 - bflying;

for bb = 1:n_tags

    landing_start    =   find(diff([blanding(1,bb);blanding(:,bb)])==1);
    landing_end      =   find(diff([blanding(:,bb);blanding(end,bb)])==-1);

    if blanding(1,bb)~=1 && blanding(end,bb)~=1
        landing_start(1) =   [];
        landing_end(1) =   [];
    elseif blanding(1,bb)==1 && blanding(end,bb)~=1
        landing_start(1) = [];
        landing_end(1:2) = [];
    elseif bflying(1,bb)~=1 && bflying(end,bb)==1
        landing_start(1) = [];
        landing_start(end) = [];
        landing_end(1) = [];
    else
        landing_start(1) = [];
        landing_start(end) = [];
        landing_end(1:2) = [];
    end

    if length(find(landing_start)) ~= length(find(landing_end))
        error('# of takeoff and landing mismatch. Check the detection.')
    end

    at_feeder{bb,1} = landing_start;
    at_feeder{bb,2} = landing_end;
    
    n_land    =   length(landing_start);
    
    feeder_distance = NaN(n_land,3);
    for ll = 1:n_land
        r_land   =   r(landing_start(ll):landing_end(ll),:,:);
        r_land_median   =   median(r_land(:,:,bb));
        
        [closest_feeder, feeder_id] = min(vecnorm(r_fd-r_land_median,2,2));
        feeder_distance(ll,1) = closest_feeder;
        feeder_distance(ll,2) = feeder_id;
        if closest_feeder < 0.5
            feeder_distance(ll,3) = 1;
        else
            feeder_distance(ll,3) = 0;
        end
    end
    at_feeder{bb,3} = feeder_distance;

    % nexttile;histogram(at_feeder{bb,3}(1,:),linspace(0,6,30))
end

%% Visualize the landing position
figure  
tiledlayout(n_tags,1,"TileSpacing","tight")
for bb = 1:n_tags
    nexttile
    plot(t,blanding(:,bb))
end

%%


    

    n_flight    =   length(takeoff_time);
    interp_pts      =   20; % Interpolate each flight into 20 points
    flight_feeder  =   zeros(n_flight,interp_pts*3+4);   % Col 1-60: interpolated position (x,y,z) during flights, col 21-24: landings on feeder 1-4

    figure
    tiledlayout(8,10,"TileSpacing","tight")
    i_figure    =   0;
    for fl  =   1:n_flight
        r_flight    =   r(takeoff_time(fl):landing_time(fl)-1,:,bb);
        r_interp    =   interp1(r_flight,linspace(1,size(r_flight,1),interp_pts),'linear');
        
        % check visually
        if i_figure < 80
            nexttile
            hold on
            plot3(r_flight(:,1),r_flight(:,2),r_flight(:,3)); view(0,90); % original
            plot3(r_interp(:,1),r_interp(:,2),r_interp(:,3)); view(0,90); % interpolated trajectory
            if i_figure == 0
                legend('original','linear interp')
            end
            hold off
    
            i_figure = i_figure + 1;
        end
        flight_feeder(fl,1:interp_pts*3)    =   reshape(r_interp,[1,interp_pts*3]);
    end
    

%% Detect feeding timings
bb   =   1;

    takeoff_time    =   find(diff([bflying(1,bb);bflying(:,bb)])==1); % time when the bat takes off
    landing_time    =   find(diff([bflying(1,bb);bflying(:,bb)])==-1); % time when the bat lands
    % Exclude the initial flight because it is the flight from a career    
    % Check if the recording does not start/end during the flight
    if bflying(1,bb)~=1 && bflying(end,bb)~=1
        takeoff_time(1) =   [];
        landing_time(1) =   [];
    elseif bflying(1,bb)==1 && bflying(end,bb)~=1
        takeoff_time(1) = [];
        landing_time(1:2) = [];
    elseif bflying(1,bb)~=1 && bflying(end,bb)==1
        takeoff_time(1) = [];
        takeoff_time(end) = [];
        landing_time(1) = [];
    else
        takeoff_time(1) = [];
        takeoff_time(end) = [];
        landing_time(1:2) = [];
    end

    if length(find(takeoff_time)) ~= length(find(landing_time))
        error('# of takeoff and landing mismatch. Check the detection.')
    end

    n_flight    =   length(takeoff_time);
    interp_pts      =   20; % Interpolate each flight into 20 points
    flight_feeder  =   zeros(n_flight,interp_pts*3+4);   % Col 1-60: interpolated position (x,y,z) during flights, col 21-24: landings on feeder 1-4

    figure
    tiledlayout(8,10,"TileSpacing","tight")
    i_figure    =   0;
    for fl  =   1:n_flight
        r_flight    =   r(takeoff_time(fl):landing_time(fl)-1,:,bb);
        r_interp    =   interp1(r_flight,linspace(1,size(r_flight,1),interp_pts),'linear');
        
        % check visually
        if i_figure < 80
            nexttile
            hold on
            plot3(r_flight(:,1),r_flight(:,2),r_flight(:,3)); view(0,90); % original
            plot3(r_interp(:,1),r_interp(:,2),r_interp(:,3)); view(0,90); % interpolated trajectory
            if i_figure == 0
                legend('original','linear interp')
            end
            hold off
    
            i_figure = i_figure + 1;
        end
        flight_feeder(fl,1:interp_pts*3)    =   reshape(r_interp,[1,interp_pts*3]);
    end

%% Correct feeder positions (s)
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position

%=== Correct feeder position(s) 
if Group_name == 'D'
    [~,fd_idx] = min(vecnorm(r_fd-centroid_all,2,2));
    r_fd(1,:) = centroid_all(fd_idx,:);                                                  
else
    for i=1:4
    [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid_all,2,2));
    if feeder_distance<0.2,r_fd(i,:) =  centroid_all(fd_idx,:);end  % Do not correct if further than 20 cm
    end
end
