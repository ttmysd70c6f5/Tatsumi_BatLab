%% Load behavioral and ephys data for all sessions
currdir = split(cd,'\');
% behdir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Behavior',currdir{end}); % behavioral data
behdir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Ephys\Ephys_Restrictive_Raw',currdir{end});  % behavioral data
ephydir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Ephys\Ephys_Restrictive_Raw',currdir{end}); % Ephys data

sessions = dir(fullfile(behdir,'**/Extracted_Behavior*.mat')); % extract paths of sessions
TTLs = dir(fullfile(ephydir,'**/Ext_Ephys_1/*c3d_1-Bat_Cluster.mat')); % reward signals

% load(fullfile(cd, 'Extracted_All.mat')) % common variables across sessions
Fs = 100;
n_tags = 6;

%% Bout detection

bmode_all = zeros(0,n_tags); % binary series of flight bouts (bouts = true) for all sessions
bout_session = zeros(0,1); % session

bout_num_all = zeros(1,n_tags); % number of bouts for all sessions
bout_len_all = cell(1,n_tags); % length of bouts
bout_sz_all = cell(1,n_tags); % number of flights in bouts
IBI_all = cell(1,n_tags); % inter bout intervals (sec)
IBI_log_all = cell(1,n_tags); % inter bout intervals (log10 sec)



% Computation of IFI
for ss = 1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying')
    [btakeoff,blanding] = findchanges(bflying,0,1);
      
    for i=1:n_tags
        % compute IFI for each session
        ifi = diff(find(btakeoff(:,i))); % IFI (sec)
        ifi_log = log10(ifi/Fs); % IFI (log10 sec)
        ifi_day = repmat(ss,[length(ifi),1]); % Experiment day
        ifi_t = find(btakeoff(:,i)); ifi_t = ifi_t(1:end-1); % timings of takeoffs
        ifi_b = repmat(i,[length(ifi),1]); % bat identification

        IFI(i).ifi       = ifi;
        IFI(i).ifi_log   = ifi_log;
        IFI(i).ifi_day   = ifi_day;
        IFI(i).ifi_t     = ifi_t;
        IFI(i).ifi_b     = ifi_b;
        
        if ss==1
            IFI_AllSession(i).ifi       = ifi;
            IFI_AllSession(i).ifi_log   = ifi_log;
            IFI_AllSession(i).ifi_day   = ifi_day;
            IFI_AllSession(i).ifi_t     = ifi_t;
            IFI_AllSession(i).ifi_b     = ifi_b;
        else
            IFI_AllSession(i).ifi       = [IFI_AllSession(i).ifi; ifi];
            IFI_AllSession(i).ifi_log   = [IFI_AllSession(i).ifi_log; ifi_log];
            IFI_AllSession(i).ifi_day   = [IFI_AllSession(i).ifi_day; ifi_day];
            IFI_AllSession(i).ifi_t     = [IFI_AllSession(i).ifi_t; ifi_t];
            IFI_AllSession(i).ifi_b     = [IFI_AllSession(i).ifi_b; ifi_b];
        end
    end
    %=== SAVE
    save(fullfile(sessions(ss).folder,'IFI.mat'),'IFI')
end
%%

%=== Bout detection and flight segmentation
load(fullfile(sessions(ss).folder,'IFI.mat'))

% % IFI clustering using k-means
% seg_method = 'k-means'; % Segmentation methods
% for i=1:n_tags
%     [ifi_idx,centroid] = kmeans(IFI_AllSession(i).ifi_log,2); % kmeans clustering
%     if centroid(1) > centroid(2) % Allocate ID=1 for smaller cluster
%         ifi_idx = 3 - ifi_idx; % Allocate cluster 1 to IFI group with smaller values
%     end
%     IFI_AllSession(i).ifi_cluster = ifi_idx;
% end

% % IFI clustering using GWM
seg_method = 'GWM'; % Segmentation methods
rng(2); % set seed
k = 2; % number of GWM clusters
option_gwm = statset('MaxIter',1000); % Options for GWM.
for i=1:n_tags    
    ifi_all = IFI_AllSession(i).ifi_log;
    gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal','SharedCovariance',false,'Options',option_gwm); % Fitted GMM
    ifi_idx = cluster(gmfit,ifi_all); % Cluster index 
    if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller cluster
        ifi_idx = 3 - ifi_idx;
    end
    IFI_AllSession(i).ifi_cluster = ifi_idx;
    
end


%=== SAVE
for ss = 1:length(sessions)
    for i=1:n_tags
        IFI(i).ifi_idx = IFI_AllSession(i).ifi_cluster(IFI_AllSession(i).ifi_day==ss);
    end
    save(fullfile(sessions(ss).folder,'IFI.mat'),'IFI')
end
%%

% Create binary timeseries data for flight bouts
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'T')

    % variables to save results for each sessions
    bout = zeros(T,n_tags);
    bout_len = cell(1,n_tags);
    bout_sz = cell(1,n_tags);
    IBI = cell(1,n_tags);
    IBI_log = cell(1,n_tags);
    
    %=== Compute flight bouts and bout sizes
    for i=1:n_tags
        ifi_now = IFI_AllSession(i).ifi(IFI_AllSession(i).ifi_day==ss);
        ifi_idx_now = IFI_AllSession(i).ifi_cluster(IFI_AllSession(i).ifi_day==ss);
        ifi_t_now = IFI_AllSession(i).ifi_t(IFI_AllSession(i).ifi_day==ss);

        bout_cnt = 0;
        bout_sz_now = []; % number of flight in a bout
        
        if ifi_t_now(find(ifi_idx_now==1,1,'first')) > 5*60*Fs % skip the first bout in initial 5 mins
            for j=1:length(ifi_now) -1
                if ifi_idx_now(j) == 1
                    bout(ifi_t_now(j):ifi_t_now(j)+ifi_now(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                    bout_cnt = bout_cnt + 1; % count of flight in bout
                    if ifi_idx_now(j+1) == 2 || j == length(ifi_now)
                        bout_sz_now = [bout_sz_now; bout_cnt];
                        bout_cnt = 1; % initialize
                    end
                end
            end
                
                
        elseif ifi_t_now(find(ifi_idx_now==1,1,'first')) <= 5*60*Fs % skip the first bout in initial 5 mins
            for j=1:length(ifi_now) -1
                if ifi_idx_now(1) == 1 % if the initial IFI belongs to the first bout, skip until the first index after the first bout
                    if j >= find(ifi_idx_now==2,1,'first') && ifi_idx_now(j) == 1
                        bout(ifi_t_now(j):ifi_t_now(j)+ifi_now(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                        bout_cnt = bout_cnt + 1; % count of flight in bout
                        if ifi_idx_now(j+1) == 2 || j == length(ifi_now)
                            bout_sz_now = [bout_sz_now; bout_cnt];
                            bout_cnt = 1; % initialize
                        end
                    end

                elseif ifi_idx_now(1) == 2 % if the initial IFI does not belong to the first bout, skip after the first bout
                    if find(ifi_idx_now==1,1,'first') < find(ifi_idx_now==2,1,'last')
                        if j >= find(ifi_idx_now(find(ifi_idx_now==1,1,'first'):end)==2,1,'first') && ifi_idx_now(j) == 1
                            bout(ifi_t_now(j):ifi_t_now(j)+ifi_now(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                            bout_cnt = bout_cnt + 1; % count of flight in bout
                            if ifi_idx_now(j+1) == 2 || j == length(ifi_now)
                                bout_sz_now = [bout_sz_now; bout_cnt];
                                bout_cnt = 1; % initialize
                            end
                        end
                    end
                end
            end
        end

        % size of bouts
        bout_day            = repmat(ss,[length(bout_sz_now),1]); % Experiment day
        [bout_takeoff,~]    = findchanges(bout,0,1);
        bout_t              = find(bout_takeoff(:,i)); % timings of initial flights of bouts
        bout_b              = repmat(i,[length(bout_sz_now),1]); % bat identification
        bout_sz{1,i} = bout_sz_now;
        bout_sz_all{1,i} = [bout_sz_all{1,i}; bout_sz_now];
    end
    
    [btakeoff,~] = findchanges(bout,0,1);
    bout_num = sum(btakeoff);
    
    bmode_all = [bmode_all; bout];
    bout_num_all = bout_num_all + bout_num;
    bout_session = [bout_session; repmat(ss,[T,1])];
    
    % Inter-bout interval
    [btakeoff,~] = findchanges(bout,0,1);
    for i=1:n_tags
         ibi = diff(find(btakeoff(:,i)))/Fs;
         ibi_log = log10(ibi);
         
         IBI{1,i} = ibi;
         IBI_log{1,i} = ibi_log;
         IBI_all{1,i} = [IBI_all{1,i}; ibi];
         IBI_log_all{1,i} = [IBI_log_all{1,i}; ibi_log];
    end
    
    % length of bouts
    for i=1:n_tags
        [btakeoff, blanding] = findchanges(bout,0,1); % timings of initial flight of bouts
        bout_len_now = (find(blanding(:,i)) - find(btakeoff(:,i)))/Fs; % length of flight bouts
        bout_len{1,i} = bout_len_now;
        bout_len_all{1,i} = [bout_len_all{1,i}; bout_len_now]; % append
    end

    
    %=== SAVE
%     save(fullfile(sessions(ss).folder,'IFI.mat'),'IFI','IFI_log','IFI_t','IFI_cluster','seg_alg')
    save(fullfile(sessions(ss).folder,'bouts.mat'),'bout','bout_num','bout_len','bout_sz','IBI','IBI_log')
    clear IFI IFI_log IFI_t IFI_cluster bout bout_num bout_len bout_sz IBI IBI_log
end

%=== SAVE
% save(fullfile(cd,'IFI_all.mat'),'ifi_AllSessions','ifi_log_AllSessions','ifi_day_AllSessions','ifi_t_AllSessions','IFI_cluster_all','IFI_combined','IFI_log_combined','seg_alg')
save(fullfile(cd,'bouts_all.mat'),'bmode_all','bout_num_all','bout_len_all','bout_sz_all','bout_session','IBI_all','IBI_log_all')
