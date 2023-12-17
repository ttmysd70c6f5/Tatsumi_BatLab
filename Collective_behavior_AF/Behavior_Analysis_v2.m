% function Behavior_Analysis_v1
% Preliminary analysis for a session
% Implement for concatenated dataset after checking with this script

%% Ephys analysis for collective behavior experiment
% Analyze d-bat data

svBhv = 'Analysis_TY/Behavior';
mkdir(svBhv)

%% Load data
%=== Ephys data
TTsFile = dir('Ext_Ephys*/SingleUnits*/*SingleUnits*');        load(fullfile(TTsFile.folder,TTsFile.name));        % TT_unit: Sorted Units
SWRFile = dir('Ext_Ephys*/CSC*/SWR*');                         load(fullfile(SWRFile.folder,SWRFile.name));        % confirmed_by_inspection: Sharp-Wave-Ripples
TTlFile = dir('Ext_Ephys*/TTL_*');                             load(fullfile(TTlFile.folder,TTlFile.name));        % TTL_timestamps: Synchronization TTLs
C3DFile = dir('Ext_Ephys*/*Bat_cluster.mat');                  load(fullfile(C3DFile.folder,C3DFile.name));        % AnalogFrameRate, AnalogSignals,Markers,VideoFrameRate: C3d File
ImpFile = dir('Ext_Ephys*/imp_bat.mat');                       load(fullfile(ImpFile.folder,ImpFile.name));        % imp_bat: number of implanted bats


%=== Behavioral data and related
BHV1_file = dir(fullfile('Ext_Beh*','Extracted_Behavior*'));                              % Preprocessed Behavioral Data
load(fullfile(BHV1_file.folder,BHV1_file.name),'a','a_abs','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
    'batdate','bflying','Fs','f_num','Group_name','n_tags','r','r_lim','t','T','v','v_abs','v_th','wBeats');        % Processed Behavioral Data

BHV2_file = dir(fullfile('Ext_Beh*','Analysis*/Analyzed_Behavior*'));
load(fullfile(BHV2_file.folder,BHV2_file.name));        % Analyzed behavior

bout_file = dir(fullfile('Ext_Beh*','bouts*'));       load(fullfile(bout_file.folder,bout_file.name));        % flight bout

bpos_file = dir(fullfile('Ext_Beh*','bpos*'));
load(fullfile(bpos_file.folder,bpos_file.name));

BHV3_file = dir(fullfile('Analysis_TY/behaviors_extracted.mat'));
load(fullfile(BHV3_file.folder,BHV3_file.name));

% Extract session number and date
dirStr = strsplit(pwd,'\');
DatasetNum = dirStr{end-1};        RecDate = dirStr{end};


%% Parameters
prm.fl_offset       =       300;                                                % Consider flights after this value (secs)

%% Summarize basic behaviors
% Use extract_behavior.m
% 
% % parameters
% fl_offset = prm.fl_offset;
% 
% % Extract flights
% bat_id          =       (1:n_tags)';
% dataset_num     =       repmat(DatasetNum,[n_tags,1]);
% rec_date        =       repmat(RecDate,[n_tags,1]);           
% fl_num          =       zeros(n_tags,1);
% 
% for nb = 1:n_tags
%     [takeoff, landing] = findchanges(bflying(:,nb),0,1);                    % Extract takeoff/landing timings
%     t_toff = t(takeoff);                                                 % Time of takeoff
%     t_land = t(landing);                                                 % Time of landing   
%     if length(t_toff) > length(t_land)
%         t_toff = t_toff(1:end-1);                                     % Remove the final takeoff in the case the session interupts the final flight
%     elseif length(t_toff) < length(t_land)
%         error('Mismatch of # of takeoffs and landings')
%     end
%     t_toff = t_toff(t_toff > fl_offset);                           % Exclude the initial flights
% 
%     fl_num(nb) = length(t_toff);                                         % Number of flights
% end
% 
% fl_rate         =       fl_num / T * Fs;                                        % Flight rate in session
% 
% % Create table and save
% behaviors = table(bat_id, dataset_num, rec_date, fl_num, fl_rate);
% save(fullfile('Analysis_TY','behaviors_extracted.mat'),'behaviors')

%% Extract all flights
% Use extract_behavior.m
% 
% t_toff      =       [];             % Time of takeoff
% t_land      =       [];             % Time of landing
% fl_len       =       [];             % Length of flight
% ts_toff     =       [];             % Timestamp of takeoff
% ts_land     =       [];             % Timestamp of landing
% bat_id      =       [];             % Bat number
% r_toff      =       [];             % location when takeoff
% r_land      =       [];             % location when landing
% pos_toff    =       [];             % position cluster when takeoff
% pos_land    =       [];         % position cluster when landing
% 
% for nb = 1:n_tags
%     [takeoff, landing] = findchanges(bflying(:,nb),0,1);                    % Extract takeoff/landing timings
%     t_toff_now       =       t(takeoff);                       
%     t_land_now       =       t(landing);                        
%     ts_toff_now      =       find(takeoff);                    
%     ts_land_now      =       find(landing);                    
%     if length(t_toff_now) > length(t_land_now)
%         t_toff_now      =       t_toff_now(1:end-1);                        % Remove the final takeoff in the case the session interupts the final flight
%         ts_toff_now     =       ts_toff_now(1:end-1);
%     elseif length(t_toff) < length(t_land)
%         error('Mismatch of # of takeoffs and landings')
%     end
% 
%     t_toff          =       [t_toff; t_toff_now];
%     t_land          =       [t_land; t_land_now];
%     fl_len           =       [fl_len; t_land_now - t_toff_now];
%     ts_toff         =       [ts_toff; ts_toff_now];
%     ts_land         =       [ts_land; ts_land_now];
%     r_toff          =       [r_toff; r(ts_toff_now,:,nb)];
%     r_land          =       [r_land; r(ts_land_now,:,nb)];
%     pos_toff        =       [pos_toff; bpos(ts_toff_now-1,nb)];
%     pos_land        =       [pos_land; bpos(ts_land_now+1,nb)];
%     bat_id          =       [bat_id; repmat(nb, [sum(takeoff),1])];
% 
% end
% 
% dataset_num     =       repmat(DatasetNum,[size(t_toff,1),1]);              % Dataset number
% rec_date        =       repmat(RecDate,[size(t_toff,1),1]);                 % Session date
% 
% % Create table and save
% flights       =       table(dataset_num, rec_date, bat_id, t_toff, t_land, ts_toff, ts_land, fl_len,...
%     r_toff, r_land, pos_toff, pos_land);
% flights       =       flights(flights.t_toff > fl_offset, :);         % Exclude initial takeoffs
% save(fullfile('Analysis_TY','behaviors_extracted.mat'),'flights','-append')

%% Extract all bouts
% Use extract_behavior.m
% 
% t_toff          =       [];             % Time of takeoff
% t_land          =       [];             % Time of landing
% bout_len          =       [];             % Length of flight
% ts_toff         =       [];             % Timestamp of takeoff
% ts_land         =       [];             % Timestamp of landing
% bat_id          =       [];             % Bat number
% r_toff          =       [];             % location when takeoff
% r_land          =       [];             % location when landing
% pos_toff        =       [];             % position cluster when takeoff
% pos_land        =       [];             % position cluster when landing
% num_fl_bout     =       [];             % Number of flights during bouts
% rate_fl_bout    =       [];             % Flight rate during bouts
% 
% for nb = 1:n_tags
%     [takeoff, landing] = findchanges(bout(:,nb),0,1);                    % Extract the first takeoff and the last landing in bouts
%     t_toff_now       =       t(takeoff);                       
%     t_land_now       =       t(landing);                        
%     ts_toff_now      =       find(takeoff);                    
%     ts_land_now      =       find(landing);                    
%     if length(t_toff_now) > length(t_land_now)
%         t_toff_now      =       t_toff_now(1:end-1);                        % Remove the final takeoff in the case the session interupts the final flight
%         ts_toff_now     =       ts_toff_now(1:end-1);
%     elseif length(t_toff) < length(t_land)
%         error('Mismatch of # of takeoffs and landings')
%     end
% 
%     [takeoff_fl, ~] = findchanges(bflying(:,nb),0,1);                   % takeoff timings of all flights
%     for nf = 1:length(t_toff_now)
%         cnt_fl       = sum(takeoff_fl(ts_toff_now(nf):ts_land_now(nf)));    % Count the number of flights during a bout
%         num_fl_bout  = [num_fl_bout; cnt_fl];
%         rate_fl_bout = [rate_fl_bout; cnt_fl ./ (t_land_now(nf)-t_toff_now(nf))];     % Compute the rate of flights during a bout
%     end
% 
%     t_toff          =       [t_toff; t_toff_now];
%     t_land          =       [t_land; t_land_now];
%     ts_toff         =       [ts_toff; ts_toff_now];
%     ts_land         =       [ts_land; ts_land_now];
%     r_toff          =       [r_toff; r(ts_toff_now,:,nb)];
%     r_land          =       [r_land; r(ts_land_now,:,nb)];
%     pos_toff        =       [pos_toff; bpos(ts_toff_now-1,nb)];
%     pos_land        =       [pos_land; bpos(ts_land_now+1,nb)];
%     bat_id          =       [bat_id; repmat(nb, [sum(takeoff),1])];
% 
% end
% 
% bout_len           =       t_land - t_toff;
% dataset_num     =       repmat(DatasetNum,[size(t_toff,1),1]);              % Dataset number
% rec_date        =       repmat(RecDate,[size(t_toff,1),1]);                 % Session date
% 
% % Create table and save
% bouts           =       table(dataset_num, rec_date, bat_id, t_toff, t_land, ts_toff, ts_land, bout_len,...
%     num_fl_bout, rate_fl_bout, r_toff, r_land, pos_toff, pos_land);
% bouts           =       bouts(bouts.t_toff > fl_offset, :);         % Exclude initial takeoffs
% save(fullfile('Analysis_TY','behaviors_extracted.mat'),'bouts','-append')

%% Spike activitiy visualization on clustered trajectories
%% PARAMETERS AND DEFINITIONS
%=== Parameters and options
n_rep = 1000;                                                                                                   % Number of repetitions for shuffling (Place Cells)
times_to_shift = t(randi([60*Fs T-60*Fs],1,n_rep));                                                           % Time intervals for circshift
n_cells = length(TT_unit);   

%=== Check if ripples were present and assign unique identifier to unit
for nc=1:n_cells
    TT_unit(nc).SWR = confirmed_by_inspection(TT_unit(nc).TT,:);
    TT_unit(nc).ID = [num2str(TT_unit(nc).bat_SN),'_', batdate,'_TT',num2str(TT_unit(nc).TT),'_', num2str(nc)];
    TT_unit(nc).An = TT_unit(nc).SWR & TT_unit(nc).Incl_flag ~=0;
    TT_unit(nc).G = Group_name;
    TT_unit(nc).Date = batdate;
end

%=== Remove cells if bad flag or no SWR
TT_unit = TT_unit([TT_unit.An]);
n_cells = length(TT_unit);

%% GENERATE THE S CELL WITH SPIKES, SURROGATE SPIKES AND AUTOCORRELATION
s = cell(n_cells,n_rep+1);      % Cell with spike times (s) between first and last TTL of the fly session + surrogates
% ac = zeros(n_cells,601);        % Autocorrelation function
Rate = zeros(length(t),n_cells);% Smoothed Firing Rate
for nc = 1:n_cells
    s{nc,1} = TT_unit(nc).Timestamps;                                                                                                     % Get the timestamps
    s{nc,1} = s{nc,1}(s{nc,1}>TTL_timestamps.fly(1) & s{nc,1}<TTL_timestamps.fly(end)); % Keep spikes between first and last session TTLs
    s{nc,1} = (s{nc,1}-TTL_timestamps.fly(1))'./1e6;  % Shift times relative to start of the fly session
    
    %=== Calculate firing rate at behavioral samples and smooth
    Rate(2:end,nc) = histcounts(s{nc,1},t)*Fs; % Histcounts return # of firings for each timepoint, so multiply Fs to get firing rate. The firing rate at time = 0 is zero.
    Rate(:,nc) = smoothdata(Rate(:,nc),1,'gaussian',round(Fs*1)); % smooth firing rate with gaussian

    for n = 1:n_rep
        s{nc,n+1} = mod(s{nc,1}+times_to_shift(n),t(T)); % Shuffling
    end        
%     [ac(nc,:),ac_xbin] = acf(s{nc,1},0.01,3);                                       % Autocorrelation
end

zRate = zscore(Rate,0,1); % z-scored firing rate

%% 3D spatial information
%% Visualization of flight trajectories and spikes
mkdir(fullfile('Analysis_TY','spike_plot'))

for nb = 1:n_tags
    for nc = 1:n_cells
        f = figure; hold on
        
        % Plot trajectories
        x = r(:,1,nb); y = r(:,2,nb);
        lh = plot(x,y,'-k');
        lh.Color(4)=0.2;
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('Spike plot: bat%d, cell%d',nb,nc))

        % Plot spikes
        sp = false(length(t),1);
        sp(2:end) = logical(histcounts(s{nc,1},t)); % 1 if the cell fires
        r_sp = r(sp,:,nb); % Position when the cell fires
        x = r_sp(:,1); y = r_sp(:,2);
        scatter(x,y,20,'r','MarkerEdgeAlpha',0.2)
        
        hold off

        exportgraphics(f,fullfile('Analysis_TY','spike_plot',sprintf('Spike_plot_bat%d_cell%d.png',nb,nc)))
        pause(0.1)
        close(gcf)
    end
end

%% Firing rate map during flights of each bat
mkdir(fullfile('Analysis_TY','rate_map'))

th_occupancy = 0.1; % minimum occupancy
bin_sz = 0.15; % spatial bin size
sigma_smooth = 1.5; % sigma for smoothing occupancy and spike counts

for nb = 1:n_tags
    f = figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.3 1 0.35]);
    tl = tiledlayout(1, n_cells, 'TileSpacing', 'tight');
    title(tl,sprintf('Firing rate map during flights: bat %d',nb))
    % nb = 6; nc = 2;

    for nc = 1:n_cells

        
        % Binning of 2D area of the room
        Xedges = r_lim(1,1)-bin_sz/3:bin_sz:r_lim(1,2)+bin_sz/3;
        Yedges = r_lim(2,1)-bin_sz/6:bin_sz:r_lim(2,2)+bin_sz/6;
        
        % Consider only flights (ignore resting)
        sp = zeros(length(t),1);
        sp(2:end) = histcounts(s{nc,1},t); % spikes: 1 if the cell fires
        sp_fly = sp(logical(sp) & logical(bflying(:,nb))); 
        r_sp = r(logical(sp) & logical(bflying(:,nb)),:,nb); % Position when the cell fires
        r_fly = r(logical(bflying(:,nb)),:,nb); % Flight trajectories
        [~,~,~,binX,binY] = histcounts2(r_sp(:,1),r_sp(:,2),Xedges,Yedges); % Histogram bin counts
        
        % Count the number of spike in each bin
        sp_cnt = zeros(length(Xedges)-1,length(Yedges)-1);
        for i = 1:size(sp_cnt,1)
            for j = 1:size(sp_cnt,2)
                sp_cnt(i,j) = sum(sp_fly(binX == i & binY ==j));
            end
        end
        % Calculate occupancy: time spent in each bin
        occupancy = histcounts2(r_fly(:,1),r_fly(:,2),Xedges,Yedges)/Fs;

        % Smooth the spike counts and the occupancy by gaussian filter
        sp_cnt = imgaussfilt(sp_cnt,sigma_smooth);
        occupancy = imgaussfilt(occupancy,sigma_smooth);

        % Invalidate bins without enough duration of stay
        occupancy_valid = false(size(occupancy));
        occupancy_pad = zeros([size(occupancy,1)+2,size(occupancy,2)+2]);
        occupancy_pad(2:end-1,2:end-1) = occupancy;
        for i = 1:size(occupancy,1)
            for j = 1:size(occupancy,2)
                occupancy_valid(i,j) = any(any(occupancy_pad(i:i+2,j:j+2) > th_occupancy));
            end
        end
        occupancy(~occupancy_valid) = nan;
        
        % Compute firing rate map by dividing the spike counts by the occupancy
        rate_map = sp_cnt ./ occupancy;
        max_fs = max(rate_map(:)); % Maximum firing rate
        % compute spatial information
        lambda = sum(sum(occupancy.*rate_map,"omitnan"));
        SI = (occupancy/sum(occupancy(:),"omitnan")).*rate_map./lambda .* log2(rate_map/lambda);
        SI = sum(SI(:),"omitnan");
        
        %===FIGURE
        nexttile
        h = imagesc(rate_map');
        set(h, 'AlphaData', 1-isnan(occupancy'))
        set(gca,'YDir','normal')
        % title(sprintf('Firing rate map: cell %d, bat %d',nc,nb))
        title(sprintf('Cell %d: %0.1fHz SI=%0.3f ',nc,max_fs,SI))
        xlabel('X')
        ylabel('Y')
        colormap jet
    end
    exportgraphics(f,fullfile('Analysis_TY','rate_map',sprintf('Firing_rate_map_bat%d.png',nb)))
    pause(0.1)
%     close gcf
end

%% Firing rate map when non-implanted bat flies but the implanted bat is resting
mkdir(fullfile('Analysis_TY','rate_map_observe'))

th_occupancy = 0.15; % minimum occupancy
sigma_smooth = 1.5; % sigma for smoothing occupancy and spike counts
bin_sz = 0.15; % spatial bin size
for nb = 1:n_tags
    if nb == imp_bat
    else
    f = figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.3 1 0.35]);
    tl = tiledlayout(1, n_cells, 'TileSpacing', 'tight');
    title(tl,sprintf('Firing rate map during observation: bat %d observed by bat %d',nb,imp_bat))
    % nb = 6; nc = 2;

    for nc = 1:n_cells

        
        % Binning of 2D area of the room
        Xedges = r_lim(1,1)-bin_sz/3:bin_sz:r_lim(1,2)+bin_sz/3;
        Yedges = r_lim(2,1)-bin_sz/6:bin_sz:r_lim(2,2)+bin_sz/6;
        
        % Consider only flights (ignore resting)
        sp = zeros(length(t),1);
        sp(2:end) = histcounts(s{nc,1},t); % spikes: 1 if the cell fires
        sp_fly = sp(logical(sp) & logical(bflying(:,nb)) & logical(1 - bflying(:,imp_bat))); 
        r_sp = r(logical(sp) & logical(bflying(:,nb)) & logical(1 - bflying(:,imp_bat)),:,nb); % Position when the cell fires
        r_fly = r(logical(bflying(:,nb) & logical(1 - bflying(:,imp_bat))),:,nb); % Flight trajectories
        [~,~,~,binX,binY] = histcounts2(r_sp(:,1),r_sp(:,2),Xedges,Yedges); % Histogram bin counts
        
        % Count the number of spike in each bin
        sp_cnt = zeros(length(Xedges)-1,length(Yedges)-1);
        for i = 1:size(sp_cnt,1)
            for j = 1:size(sp_cnt,2)
                sp_cnt(i,j) = sum(sp_fly(binX == i & binY ==j));
            end
        end
        % Calculate occupancy: time spent in each bin
        occupancy = histcounts2(r_fly(:,1),r_fly(:,2),Xedges,Yedges)/Fs;
        
        % Smooth the spike counts and the occupancy by gaussian filter
        sp_cnt = imgaussfilt(sp_cnt,sigma_smooth);
        occupancy = imgaussfilt(occupancy,sigma_smooth);

        % Invalidate bins without enough duration of stay
        occupancy_valid = false(size(occupancy));
        occupancy_pad = zeros([size(occupancy,1)+2,size(occupancy,2)+2]);
        occupancy_pad(2:end-1,2:end-1) = occupancy;
        for i = 1:size(occupancy,1)
            for j = 1:size(occupancy,2)
                occupancy_valid(i,j) = any(any(occupancy_pad(i:i+2,j:j+2) > th_occupancy));
            end
        end
        occupancy(~occupancy_valid) = nan;
        
        % Compute firing rate map by dividing the spike counts by the occupancy
        rate_map = sp_cnt ./ occupancy;
        max_fs = max(rate_map(:)); % Maximum firing rate
        % compute spatial information
        lambda = sum(sum(occupancy.*rate_map,"omitnan"));
        SI = (occupancy/sum(occupancy(:),"omitnan")).*rate_map./lambda .* log2(rate_map/lambda);
        SI = sum(SI(:),"omitnan");
        
        %===FIGURE
        nexttile
        h = imagesc(rate_map');
        set(h, 'AlphaData', 1-isnan(occupancy'))
        set(gca,'YDir','normal')
        % title(sprintf('Firing rate map: cell %d, bat %d',nc,nb))
        title(sprintf('Cell %d: %0.1fHz SI=%0.3f ',nc,max_fs,SI))
        xlabel('X')
        ylabel('Y')
        colormap jet
    end
    exportgraphics(f,fullfile('Analysis_TY','rate_map_observe',sprintf('Firing_rate_map_bat%d.png',nb)))
    pause(0.1)
%     close gcf
    end
end


%% Firing rate map when implanted bat flies but non-implanted bat is resting
mkdir(fullfile('Analysis_TY','rate_map_observe'))

th_occupancy = 0.15; % minimum occupancy
sigma_smooth = 1.5; % sigma for smoothing occupancy and spike counts
bin_sz = 0.15; % spatial bin size

nb = imp_bat;

f = figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.3 1 0.35]);
tl = tiledlayout(1, n_cells, 'TileSpacing', 'tight');
title(tl,sprintf('Firing rate map: bat %d is flying and other bats are resting',imp_bat))

for nc = 1:n_cells
    % Binning of 2D area of the room
    Xedges = r_lim(1,1)-bin_sz/3:bin_sz:r_lim(1,2)+bin_sz/3;
    Yedges = r_lim(2,1)-bin_sz/6:bin_sz:r_lim(2,2)+bin_sz/6;
    
    % Consider only flights (ignore resting)
    sp = zeros(length(t),1);
    sp(2:end) = histcounts(s{nc,1},t); % spikes: 1 if the cell fires
    sp_fly = sp(logical(sp) & any(logical(1-bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat))); 
    r_sp = r(logical(sp) & any(logical(1-bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat)),:,nb); % Position when the cell fires
    r_fly = r(any(logical(1-bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat)),:,nb); % Flight trajectories
    [~,~,~,binX,binY] = histcounts2(r_sp(:,1),r_sp(:,2),Xedges,Yedges); % Histogram bin counts
    
    % Count the number of spike in each bin
    sp_cnt = zeros(length(Xedges)-1,length(Yedges)-1);
    for i = 1:size(sp_cnt,1)
        for j = 1:size(sp_cnt,2)
            sp_cnt(i,j) = sum(sp_fly(binX == i & binY ==j));
        end
    end
    % Calculate occupancy: time spent in each bin
    occupancy = histcounts2(r_fly(:,1),r_fly(:,2),Xedges,Yedges)/Fs;
    
    % Smooth the spike counts and the occupancy by gaussian filter
    sp_cnt = imgaussfilt(sp_cnt,sigma_smooth);
    occupancy = imgaussfilt(occupancy,sigma_smooth);

    % Invalidate bins without enough duration of stay
    occupancy_valid = false(size(occupancy));
    occupancy_pad = zeros([size(occupancy,1)+2,size(occupancy,2)+2]);
    occupancy_pad(2:end-1,2:end-1) = occupancy;
    for i = 1:size(occupancy,1)
        for j = 1:size(occupancy,2)
            occupancy_valid(i,j) = any(any(occupancy_pad(i:i+2,j:j+2) > th_occupancy));
        end
    end
    occupancy(~occupancy_valid) = nan;
    
    % Compute firing rate map by dividing the spike counts by the occupancy
    rate_map = sp_cnt ./ occupancy;
    max_fs = max(rate_map(:)); % Maximum firing rate
    % compute spatial information
    lambda = sum(sum(occupancy.*rate_map,"omitnan"));
    SI = (occupancy/sum(occupancy(:),"omitnan")).*rate_map./lambda .* log2(rate_map/lambda);
    SI = sum(SI(:),"omitnan");
    
    %===FIGURE
    nexttile
    h = imagesc(rate_map');
    set(h, 'AlphaData', 1-isnan(occupancy'))
    set(gca,'YDir','normal')
    % title(sprintf('Firing rate map: cell %d, bat %d',nc,nb))
    title(sprintf('Cell %d: %0.1fHz SI=%0.3f ',nc,max_fs,SI))
    xlabel('X')
    ylabel('Y')
    colormap jet
end
exportgraphics(f,fullfile('Analysis_TY','rate_map_observe',sprintf('Firing_rate_map_self.png')))
pause(0.1)
close gcf

%% Firing rate map when both implanted bat and non-implanted bats fly
% th_occupancy = 0.1; % minimum occupancy
% sigma_smooth = 1.5; % sigma for smoothing occupancy and spike counts
% bin_sz = 0.15; % spatial bin size
% 
% nb = imp_bat;
% 
% f = figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.3 1 0.35]);
% tl = tiledlayout(1, n_cells, 'TileSpacing', 'tight');
% title(tl,sprintf('Firing rate map: bat %d is flying and other bats are resting',imp_bat))
% 
% for nc = 1:n_cells
%     % Binning of 2D area of the room
%     Xedges = r_lim(1,1)-bin_sz/3:bin_sz:r_lim(1,2)+bin_sz/3;
%     Yedges = r_lim(2,1)-bin_sz/6:bin_sz:r_lim(2,2)+bin_sz/6;
%     
%     % Consider only flights (ignore resting)
%     sp = zeros(length(t),1);
%     sp(2:end) = histcounts(s{nc,1},t); % spikes: 1 if the cell fires
%     sp_fly = sp(logical(sp) & any(logical(bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat))); 
%     r_sp = r(logical(sp) & any(logical(bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat)),:,nb); % Position when the cell fires
%     r_fly = r(any(logical(bflying(:,[1:imp_bat-1,imp_bat+1:n_tags])),2) & logical(bflying(:,imp_bat)),:,nb); % Flight trajectories
%     [~,~,~,binX,binY] = histcounts2(r_sp(:,1),r_sp(:,2),Xedges,Yedges); % Histogram bin counts
%     
%     % Count the number of spike in each bin
%     sp_cnt = zeros(length(Xedges)-1,length(Yedges)-1);
%     for i = 1:size(sp_cnt,1)
%         for j = 1:size(sp_cnt,2)
%             sp_cnt(i,j) = sum(sp_fly(binX == i & binY ==j));
%         end
%     end
%     % Calculate occupancy: time spent in each bin
%     occupancy = histcounts2(r_fly(:,1),r_fly(:,2),Xedges,Yedges)/Fs;
%     
%     % Smooth the spike counts and the occupancy by gaussian filter
%     sp_cnt = imgaussfilt(sp_cnt,sigma_smooth);
%     occupancy = imgaussfilt(occupancy,sigma_smooth);
% 
%     % Invalidate bins without enough duration of stay
%     occupancy_valid = false(size(occupancy));
%     occupancy_pad = zeros([size(occupancy,1)+2,size(occupancy,2)+2]);
%     occupancy_pad(2:end-1,2:end-1) = occupancy;
%     for i = 1:size(occupancy,1)
%         for j = 1:size(occupancy,2)
%             occupancy_valid(i,j) = any(any(occupancy_pad(i:i+2,j:j+2) > th_occupancy));
%         end
%     end
%     occupancy(~occupancy_valid) = nan;
%     
%     % Compute firing rate map by dividing the spike counts by the occupancy
%     rate_map = sp_cnt ./ occupancy;
%     max_fs = max(rate_map(:)); % Maximum firing rate
%     % compute spatial information
%     lambda = sum(sum(occupancy.*rate_map,"omitnan"));
%     SI = (occupancy/sum(occupancy(:),"omitnan")).*rate_map./lambda .* log2(rate_map/lambda);
%     SI = sum(SI(:),"omitnan");
%     
%     %===FIGURE
%     nexttile
%     h = imagesc(rate_map');
%     set(h, 'AlphaData', 1-isnan(occupancy'))
%     set(gca,'YDir','normal')
%     % title(sprintf('Firing rate map: cell %d, bat %d',nc,nb))
%     title(sprintf('Cell %d: %0.1fHz SI=%0.3f ',nc,max_fs,SI))
%     xlabel('X')
%     ylabel('Y')
%     colormap jet
% end

%% Trajectory clustering
n_smp = 20; % number of downsampled points
n_dim = 3; % number of dimension (x,y,z)
cutoff = 1.1; % Cutoff value for clustering
min_trj = 3; % Minimum amounts of trajectory in a cluster

r_interp = zeros(size(flights,1),n_dim*n_smp); % interpolated trajectory
sp_interp = zeros(size(flights,1),n_smp,n_cells);
for nf = 1:size(flights,1)
    nb = flights.bat_id(nf); % Bat id
    r_flight = r(flights.ts_toff(nf):flights.ts_land(nf),:,nb); % raw trajectory
    sp_flight = Rate(flights.ts_toff(nf):flights.ts_land(nf),:); % firing rate
    fl_len_ts = flights.ts_land(nf) - flights.ts_toff(nf); % length of flight (timestamps)
    w_interp = 1 + fl_len_ts .* (0.5 + 0.5 .* (linspace(-n_smp/2,n_smp/2,n_smp)).^3 / (n_smp/2)^3); % weight for interpolation: sigmoid
%     w_interp = linspace(1,size(r_flight,1),n_smp);
    r_interp(nf,:) = reshape(interp1(r_flight, w_interp),[1,n_smp*3]);
    sp_interp(nf,:,:) = interp1(sp_flight,w_interp);
end

%% Hierarchical clustering and visualize rate map
mkdir('Analysis_TY/hc')
hc = struct([]);
min_trj = 5;
nb = 6;
% hierarchical clustering
hc_tree = linkage(r_interp(flights.bat_id==nb,:),'single','euclidean'); % hierarchical cluster tree
D = pdist(r_interp(flights.bat_id==nb,:));
leafOrder = optimalleaforder(hc_tree,D);
hc(nb).hc_clus = cluster(hc_tree,'cutoff',cutoff,'Criterion','inconsistent');
[hc(nb).hc_cts,hc(nb).hc_id] = groupcounts(hc(nb).hc_clus);
hc(nb).hc_good = hc(nb).hc_id(hc(nb).hc_cts > min_trj);

fprintf('\n Number of clusters for bat %d: %d \n',nb,length(unique(hc(nb).hc_id)))
fprintf('Number of clusters with more than %d flights for bat %d: %d \n',min_trj,nb,length(unique(hc(nb).hc_good)))

% takeoff and landing position for good clusters
for nclu = 1:length(hc(nb).hc_good)
    r_toff_ext = flights.r_toff(flights.bat_id==nb,:);
    r_land_ext = flights.r_land(flights.bat_id==nb,:);

    hc(nb).r_toff_mean(nclu,:) = mean(r_toff_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
    hc(nb).r_land_mean(nclu,:) = mean(r_land_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
end

% sort good clusters based on its takeoff and landing positions
toff_land_pos = [hc(nb).r_toff_mean, hc(nb).r_land_mean];
toff_land_pos = round(toff_land_pos);
[~,sorted_idx] = sortrows(toff_land_pos,[1,2,4,5]);
hc(nb).hc_good = hc(nb).hc_good(sorted_idx);
hc(nb).r_toff_mean = hc(nb).r_toff_mean(sorted_idx);
hc(nb).r_land_mean = hc(nb).r_land_mean(sorted_idx);

%===FIGURE: visualize clustered trajectories for each bat
f = figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl = tiledlayout(3,5,'TileSpacing','tight');

flights_ext = flights(flights.bat_id==nb,:); % flights
r_interp_ext = r_interp(flights.bat_id==nb,:); % position
r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);

for nclu = 1:length(hc(nb).hc_good)
    nexttile; hold on

    flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
    r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:,:);
    
%         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
%         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
    c = 1:.1:n_smp; c = c/max(c); c = 100.^c;
    hold on
    for nf = 1:size(flights_clu_now,1)
        x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
        y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
        scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
    end
    hold off
    colormap parula
    xlim(r_lim(1,:)); ylim(r_lim(2,:));
    title(tl,sprintf('Hierarchical clustering of trajectory: cell %d, bat %d, cutoff=%0.1f, minimum number of trajectory=%d',nc,nb,round(cutoff,1),min_trj))
    xlabel(tl,'x'); ylabel(tl,'y')
end
exportgraphics(f,fullfile('Analysis_TY','hc',sprintf('Trajectory_hc_bat%d.png',nb)))
pause(0.1)
% close gcf

%% 
nclu = 2;
flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
r_clu_now =  flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
nc = 1;
bin_sz = 1.5;

% Binning of 2D area of the room
Xedges = r_lim(1,1)-bin_sz/3:bin_sz:r_lim(1,2)+bin_sz/3;
Yedges = r_lim(2,1)-bin_sz/6:bin_sz:r_lim(2,2)+bin_sz/6;

% Consider only flights (ignore resting)
sp = zeros(length(t),1);
sp(2:end) = histcounts(s{nc,1},t); % spikes: 1 if the cell fires
sp_fly = sp(logical(sp) & logical(bflying(:,nb))); 
r_sp = r(logical(sp) & logical(bflying(:,nb)),:,nb); % Position when the cell fires
r_fly = r(logical(bflying(:,nb)),:,nb); % Flight trajectories
[~,~,~,binX,binY] = histcounts2(r_sp(:,1),r_sp(:,2),Xedges,Yedges); % Histogram bin counts

% Count the number of spike in each bin
sp_cnt = zeros(length(Xedges)-1,length(Yedges)-1);
for i = 1:size(sp_cnt,1)
    for j = 1:size(sp_cnt,2)
        sp_cnt(i,j) = sum(sp_fly(binX == i & binY ==j));
    end
end
% Calculate occupancy: time spent in each bin
occupancy = histcounts2(r_fly(:,1),r_fly(:,2),Xedges,Yedges)/Fs;
occupancy_valid = false(size(occupancy));
occupancy_pad = zeros([size(occupancy,1)+2,size(occupancy,2)+2]);
occupancy_pad(2:end-1,2:end-1) = occupancy;
for i = 1:size(occupancy,1)
    for j = 1:size(occupancy,2)
        occupancy_valid(i,j) = any(any(occupancy_pad(i:i+2,j:j+2) > th_occupancy));
    end
end
occupancy(~occupancy_valid) = nan;
% Smooth the spike counts and the occupancy by gaussian filter
sp_cnt = imgaussfilt(sp_cnt,1);
occupancy = imgaussfilt(occupancy,1);
% Compute firing rate map by dividing the spike counts by the occupancy
rate_map = sp_cnt ./ occupancy;
max_fs = max(rate_map(:)); % Maximum firing rate
% compute spatial information
lambda = sum(sum(occupancy.*rate_map,"omitnan"));
SI = (occupancy/sum(occupancy(:),"omitnan")).*rate_map./lambda .* log2(rate_map/lambda);
SI = sum(SI(:),"omitnan");

%===FIGURE
nexttile
h = imagesc(rate_map');
set(h, 'AlphaData', 1-isnan(occupancy'))
set(gca,'YDir','normal')
% title(sprintf('Firing rate map: cell %d, bat %d',nc,nb))
title(sprintf('Cell %d: %0.1fHz SI=%0.3f ',nc,max_fs,SI))
xlabel('X')
ylabel('Y')
colormap jet

%% hierarchical clustering
hc = struct([]);

for nb=1:n_tags
    % hierarchical clustering
    hc_tree = linkage(r_interp(flights.bat_id==nb,:),'single','euclidean'); % hierarchical cluster tree
    D = pdist(r_interp(flights.bat_id==nb,:));
    leafOrder = optimalleaforder(hc_tree,D);
    hc(nb).hc_clus = cluster(hc_tree,'cutoff',cutoff,'Criterion','inconsistent');
    [hc(nb).hc_cts,hc(nb).hc_id] = groupcounts(hc(nb).hc_clus);
    hc(nb).hc_good = hc(nb).hc_id(hc(nb).hc_cts > min_trj);

    fprintf('\n Number of clusters for bat %d: %d \n',nb,length(unique(hc(nb).hc_id)))
    fprintf('Number of clusters with more than %d flights for bat %d: %d \n',min_trj,nb,length(unique(hc(nb).hc_good)))

    % takeoff and landing position for good clusters
    for nclu = 1:length(hc(nb).hc_good)
        r_toff_ext = flights.r_toff(flights.bat_id==nb,:);
        r_land_ext = flights.r_land(flights.bat_id==nb,:);

        hc(nb).r_toff_mean(nclu,:) = mean(r_toff_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
        hc(nb).r_land_mean(nclu,:) = mean(r_land_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
    end

    % sort good clusters based on its takeoff and landing positions
    toff_land_pos = [hc(nb).r_toff_mean, hc(nb).r_land_mean];
    toff_land_pos = round(toff_land_pos);
    [~,sorted_idx] = sortrows(toff_land_pos,[1,2,4,5]);
    hc(nb).hc_good = hc(nb).hc_good(sorted_idx);
    hc(nb).r_toff_mean = hc(nb).r_toff_mean(sorted_idx);
    hc(nb).r_land_mean = hc(nb).r_land_mean(sorted_idx);

    %===FIGURE: visualize clustered trajectories for each bat
    f = figure;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    tl = tiledlayout(3,5,'TileSpacing','tight');

    flights_ext = flights(flights.bat_id==nb,:); % flights
    r_interp_ext = r_interp(flights.bat_id==nb,:); % position
    r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);

    for nclu = 1:length(hc(nb).hc_good)
        nexttile; hold on

        flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
        r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:,:);
        
%         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
%         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
        c = 1:.1:n_smp; c = c/max(c); c = 100.^c;
        hold on
        for nf = 1:size(flights_clu_now,1)
            x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
            y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
            scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
        end

        hold off
        colormap parula
        xlim(r_lim(1,:)); ylim(r_lim(2,:));
        title(tl,sprintf('Hierarchical clustering of trajectory: cell %d, bat %d, cutoff=%0.1f, minimum number of trajectory=%d',nc,nb,round(cutoff,1),min_trj))
        xlabel(tl,'x'); ylabel(tl,'y')
    end
    exportgraphics(f,fullfile('Analysis_TY','hc',sprintf('Trajectory_hc_bat%d.png',nb)))
    pause(0.1)
    close gcf



    %===FIGURE: firing rate on trajectories
    for nc = 1:n_cells
        f = figure;
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        tl = tiledlayout(3,5,'TileSpacing','tight');
    
        flights_ext = flights(flights.bat_id==nb,:);
        r_interp_ext = r_interp(flights.bat_id==nb,:);
        r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);
        sp_interp_ext    = sp_interp(flights.bat_id==nb,:,nc); % firing rate

        for nclu = 1:length(hc(nb).hc_good)
            nexttile; hold on
    
            flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
            r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:,:);
            sp_interp_clu_now = sp_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
            
    %         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
    %         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
        
            for nf = 1:size(flights_clu_now,1)
                x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
                y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
%                 c = spline(1:n_smp,sp_interp_clu_now(nf,:),1:.1:n_smp); 
%                 if max(c) ~= 0
%                     c = c/max(c); c = 100.^c;
%                 end              
%                 scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)  
                plot(x,y,'-k')
            end

            % Plot spikes
            sp = false(length(t),1);
            sp(2:end) = logical(histcounts(s{nc,1},t)); % 1 if the cell fires
            r_sp = r(sp,:,nb); % Position when the cell fires
        x = r_sp(:,1); y = r_sp(:,2);
        scatter(x,y,20,'r','MarkerEdgeAlpha',0.2)

            hold off
            colormap jet
            xlim(r_lim(1,:)); ylim(r_lim(2,:));
            title(tl,sprintf('Firing rate on trajectory: cell %d, bat %d, cutoff=%0.1f, minimum number of trajectory=%d',nc,nb,round(cutoff,1),min_trj))
            xlabel(tl,'x'); ylabel(tl,'y')
        end
        exportgraphics(f,fullfile('Analysis_TY','hc',sprintf('Firing_hc_bat%d_cell%d.png',nb,nc)))
        
        pause(0.1)
        close gcf
    end
end

%% Analysis for distance from other bat
r_dist = zeros(size(r,1),size(r,2),n_tags-1);
for nb = 1:n_tags-1
    r_dist(:,:,nb) = r(:,:,nb) - r(:,:,imp_bat);
end

% Distance clustering
n_smp = 20; % number of downsampled points
n_dim = 3; % number of dimension (x,y,z)
cutoff = 1.1; % Cutoff value for clustering
min_trj = 3; % Minimum amounts of trajectory in a cluster

flights_noimp = flights(flights.bat_id~=imp_bat,:);
r_interp = zeros(size(flights_noimp,1),n_dim*n_smp); % interpolated trajectory
sp_interp = zeros(size(flights_noimp,1),n_smp,n_cells);
for nf = 1:size(flights_noimp,1)
    nb = flights_noimp.bat_id(nf); % Bat id
    r_flight = r_dist(flights_noimp.ts_toff(nf):flights_noimp.ts_land(nf),:,nb); % raw trajectory
    sp_flight = Rate(flights_noimp.ts_toff(nf):flights_noimp.ts_land(nf),:); % firing rate
    fl_len_ts = flights_noimp.ts_land(nf) - flights_noimp.ts_toff(nf); % length of flight (timestamps)
    w_interp = 1 + fl_len_ts .* (0.5 + 0.5 .* (linspace(-n_smp/2,n_smp/2,n_smp)).^3 / (n_smp/2)^3); % weight for interpolation: sigmoid
%     w_interp = linspace(1,size(r_flight,1),n_smp);
    r_interp(nf,:) = reshape(interp1(r_flight, w_interp),[1,n_smp*3]);
    sp_interp(nf,:,:) = interp1(sp_flight,w_interp);
end

hc = struct([]);

for nb=1:n_tags-1
    % hierarchical clustering
    hc_tree = linkage(r_interp(flights_noimp.bat_id==nb,:),'single','euclidean'); % hierarchical cluster tree
    D = pdist(r_interp(flights_noimp.bat_id==nb,:));
    leafOrder = optimalleaforder(hc_tree,D);
    hc(nb).hc_clus = cluster(hc_tree,'cutoff',cutoff,'Criterion','inconsistent');
    [hc(nb).hc_cts,hc(nb).hc_id] = groupcounts(hc(nb).hc_clus);
    hc(nb).hc_good = hc(nb).hc_id(hc(nb).hc_cts > min_trj);

    fprintf('\n Number of clusters for bat %d: %d \n',nb,length(unique(hc(nb).hc_id)))
    fprintf('Number of clusters with more than %d flights for bat %d: %d \n',min_trj,nb,length(unique(hc(nb).hc_good)))

    % takeoff and landing position for good clusters
    for nclu = 1:length(hc(nb).hc_good)
        r_toff_ext = flights_noimp.r_toff(flights_noimp.bat_id==nb,:);
        r_land_ext = flights_noimp.r_land(flights_noimp.bat_id==nb,:);

        hc(nb).r_toff_mean(nclu,:) = mean(r_toff_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
        hc(nb).r_land_mean(nclu,:) = mean(r_land_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:));
    end

    % sort good clusters based on its takeoff and landing positions
    toff_land_pos = [hc(nb).r_toff_mean, hc(nb).r_land_mean];
    toff_land_pos = round(toff_land_pos);
    [~,sorted_idx] = sortrows(toff_land_pos,[1,2,4,5]);
    hc(nb).hc_good = hc(nb).hc_good(sorted_idx);
    hc(nb).r_toff_mean = hc(nb).r_toff_mean(sorted_idx);
    hc(nb).r_land_mean = hc(nb).r_land_mean(sorted_idx);

    %===FIGURE: visualize clustered trajectories for each bat
    f = figure;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    tl = tiledlayout(3,5,'TileSpacing','tight');

    flights_ext = flights_noimp(flights_noimp.bat_id==nb,:); % flights
    r_interp_ext = r_interp(flights_noimp.bat_id==nb,:); % position
    r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);

    for nclu = 1:length(hc(nb).hc_good)
        nexttile; hold on

        flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
        r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:,:);
        
%         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
%         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
        c = 1:.1:n_smp; c = c/max(c); c = 100.^c;
        hold on
        for nf = 1:size(flights_clu_now,1)
            x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
            y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
            scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
        end
        hold off
        colormap parula
        xlim(r_lim(1,:)); ylim(r_lim(2,:));
        title(tl,sprintf('Hierarchical clustering of distance: cell %d, bat %d vs imp bat, cutoff=%0.1f, minimum number of trajectory=%d',nc,nb,round(cutoff,1),min_trj))
        xlabel(tl,'x'); ylabel(tl,'y')
    end
    exportgraphics(f,fullfile('Analysis_TY',sprintf('Dist_hc_bat%d_imp.png',nb)))
    pause(0.1)
    close gcf

    %===FIGURE: firing rate on trajectories
    for nc = 1:n_cells
        f = figure;
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        tl = tiledlayout(3,5,'TileSpacing','tight');
    
        flights_ext = flights_noimp(flights_noimp.bat_id==nb,:);
        r_interp_ext = r_interp(flights_noimp.bat_id==nb,:);
        r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);
        sp_interp_ext    = sp_interp(flights_noimp.bat_id==nb,:,nc); % firing rate

        for nclu = 1:length(hc(nb).hc_good)
            nexttile; hold on
    
            flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
            r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:,:);
            sp_interp_clu_now = sp_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nclu),:);
            
    %         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
    %         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
        
            for nf = 1:size(flights_clu_now,1)
                x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
                y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
                c = spline(1:n_smp,sp_interp_clu_now(nf,:),1:.1:n_smp); 
                if max(c) ~= 0
                    c = c/max(c); c = 100.^c;
                end              
                scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
            end
            hold off
            colormap jet
            xlim(r_lim(1,:)); ylim(r_lim(2,:));
            title(tl,sprintf('Firing rate vs distance: bat %d vs imp, cutoff=%0.1f, minimum number of trajectory=%d',nb,round(cutoff,1),min_trj))
            xlabel(tl,'x'); ylabel(tl,'y')
        end
        exportgraphics(f,fullfile('Analysis_TY',sprintf('Dist_firing_hc_bat%d_imp_cell%d.png',nb,nc)))
        
        pause(0.1)
        close gcf
    end
end

