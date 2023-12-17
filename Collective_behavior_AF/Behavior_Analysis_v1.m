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
n_rep = 100;                                                                                                   % Number of repetitions for shuffling (Place Cells)
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
    exportgraphics(f,fullfile('Analysis_TY',sprintf('Trajectory_hc_bat%d.png',nb)))
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
                c = spline(1:n_smp,sp_interp_clu_now(nf,:),1:.1:n_smp); 
                if max(c) ~= 0
                    c = c/max(c); c = 100.^c;
                end              
                scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
            end
            hold off
            colormap jet
            xlim(r_lim(1,:)); ylim(r_lim(2,:));
            title(tl,sprintf('Firing rate on trajectory: bat %d, cutoff=%0.1f, minimum number of trajectory=%d',nb,round(cutoff,1),min_trj))
            xlabel(tl,'x'); ylabel(tl,'y')
        end
        exportgraphics(f,fullfile('Analysis_TY',sprintf('Firing_hc_bat%d_cell%d.png',nb,nc)))
        
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

