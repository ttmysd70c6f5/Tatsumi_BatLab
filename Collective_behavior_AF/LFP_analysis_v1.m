%% Function for aggregating and filtering raw voltage traces from the same Tetrode into LFP
% The function runs for each tetrode by averaging the corresponding channels, aligning in time,
% filtering between 0.5 and 300 Hz and downsampling at 1 kHz. 

%=== Load Collateral Files
load('TTL_timestamps.mat');                                                             % Reference Timestamps for alignment
load('imp_bat.mat');                                                                    % Tag of the implanted bat
BHV_file = dir(fullfile(fileparts(cd),'Ext_Beh*','Extracted_Behavior*'));               % Behavioral variables (velocity, acceleration, etc...)
load([BHV_file.folder,'\',BHV_file.name],'r','v_abs','a_abs','t','bflying','batdate','n_tags','T');

%=== Get the name of the folder for saving the data
analysis_directory=fullfile(pwd,['CSC_Analysis_',batdate,]);

%=== Options and Init
options.savefigures = 1;        % Save Figures
fig_count = 1;                  % Id of the first figure
TT_LFP_channels = zeros(4);     % Save excluded tetrode channels
Freq_band = [0.5 300];          % Frequency band for filtering the signal

%=== Process channels for each of the 4 Tetrodes
for TT = 1:4
    
    %=== Define channels corresponding to the TT and initialize
    disp(['Analyzing voltage from TT',num2str(TT), ' ...']);
    CSC = 4*(TT-1)+[0:3];
    
    for ch = 1:4
        %=== Load CSC Files and assign to variables
        CSC_file = dir(fullfile(cd, ['*CSC' num2str(CSC(ch)),'.mat']));
        load(CSC_file.name,'Estimated_channelFS_Transceiver','AD_count_int16','AD_count_to_uV_factor','Timestamps_of_first_samples_usec');
        Fs_CSC = Estimated_channelFS_Transceiver(1);                                                % Voltage Sampling Frequency (Hz)
        raw_V{ch,:} = double(AD_count_int16*AD_count_to_uV_factor);                                 % Raw Voltage Trace
        tstamps{ch,:} = Timestamps_of_first_samples_usec(1)+[0:length(raw_V{ch,:})-1]/Fs_CSC*1e6;   % Timestamps (us)
    end
    
    %=== Calculate differences to spot weird channels
    channel_pairs = nchoosek(1:4,2);
    for i=1:size(channel_pairs,1)
    channel_pairs(i,3) = std(abs(raw_V{channel_pairs(i,1),:}-raw_V{channel_pairs(i,2),:}));
    end
    
    %=== Assume there is one bad channel
    bad_flag = logical(kmeans(channel_pairs(:,3),2)-1);    % assume 0 are the low difference pairs
    if mean(channel_pairs(~bad_flag,3))>mean(channel_pairs(bad_flag,3))
        bad_flag = ~bad_flag;
    end
   
    %=== Remove the channel if the model is consistent (3 pairs containing the same channel are bad)
    good_channels = true(4,1);
    if mean(channel_pairs(bad_flag,3))/mean(channel_pairs(~bad_flag,3))>5 && nnz(bad_flag)==3 
        i = 1;
        while ~all(any(channel_pairs(bad_flag,1:2) == i,2)) && i<=4
        i=i+1;
        end
        disp(['Removing Ch',num2str(i)]);
        good_channels(i) = 0;
    end
    
    %=== Calculate average of good channels
    mean_V = mean(vertcat(raw_V{good_channels,1}),1);
    TT_LFP_channels(TT,:) = good_channels;
    
    %=== Check timestamps
    ts_differences = zeros(1,6);
    for i=1:size(channel_pairs,1)
        ts_differences(i) = mean(abs(tstamps{channel_pairs(i,1),:}-tstamps{channel_pairs(i,2),:}));
    end
    if ~all(ts_differences<1e3)
        disp('--------------------Large timestamp differences detected------------------!!!!!!');
    else
        disp('Timestamps are good');
    end
    mean_tstamps = tstamps{1,:};
    
    %=== Perform temporal alignment
    [~,sample_strt] = min(abs(mean_tstamps-TTL_timestamps.fly(1)));                     % First timestamp of the fly session
    [~,sample_stop] = min(abs(mean_tstamps-TTL_timestamps.fly(end)));                   % Last timestamp of the fly session
    V_fly  = mean_V(:,sample_strt:sample_stop);                                         % Voltage (averaged across good channels)
    t_fly = mean_tstamps(:,sample_strt:sample_stop);   t_fly = (t_fly-t_fly(1))/1e6;    % Time in seconds
    
    %=== Downsample time to 1 Khz (no need to keep higher sampling)
    t_1KHz = t_fly(1):1/1e3:t_fly(end);
    
    %=== Low pass filter between 0.5 and 300 Hz and then interpolate at downsampled time
    disp('Filtering LFP...');
    LFP_temp = [];
    LFP_temp = interp1(t_fly,bandpass(V_fly,Freq_band,Fs_CSC,ImpulseResponse="iir"),t_1KHz);   % 0.85 steepness, 60 dB
    
    %=== Adjust for rare differences in size across tetrodes
    if TT>1 && numel(LFP_temp)>numel(LFP(1,:))
        LFP_temp = LFP_temp(1:numel(LFP(1,:)));
    elseif TT>1 && numel(LFP_temp)<numel(LFP(1,:))
        LFP_temp = [LFP_temp,LFP_temp(end)*ones(1,numel(LFP(1,:))-numel(LFP_temp))];
    end
    
    %Assign to the LFP matrix
    LFP(TT,:) = LFP_temp;   

    %=== Show a sample of the session
    t1 = 150;   t2 = t1+10;
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(5,2,'TileSpacing','tight');
    for i=1:10
        nexttile;   plot(t_1KHz(t_1KHz>i*t1 & t_1KHz<i*t1+10),LFP(TT,t_1KHz>i*t1 & t_1KHz<i*t1+10),'r');   hold on;
                    plot(t_fly(t_fly>i*t1 & t_fly<i*t1+10),V_fly(t_fly>i*t1 & t_fly<i*t1+10),'k'); ylim([-1e3 1e3]);
    xlabel('Time (s)'); ylabel('Voltage (uV)');
    end
    sgtitle(['LFP, TT', num2str(TT), ' , ', num2str(batdate)]);
    fig_count = saveFig(analysis_directory,['LFP_', batdate],fig_count,options.savefigures);
    
    %plot(raw_V{1,:}(5e4:6e4));    hold on;  plot(raw_V{2,:}(5e4:6e4));    plot(raw_V{3,:}(5e4:6e4));     plot(raw_V{4,:}(5e4:6e4)); plot(mean_V(5e4:6e4),'k');   
  
end

%=== Save the data
LFP_data.LFP = LFP;
LFP_data.t_1KHz = t_1KHz;
LFP_data.TT_LFP_channels = TT_LFP_channels;
LFP_data.Freq_band = Freq_band;
save([analysis_directory,'/LFP_Ext_', batdate,'.mat'],'LFP_data');