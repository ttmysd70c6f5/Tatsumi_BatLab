% function Ephys_Analysis_v1(ephysdir)
% currdir = pwd;
% cd(ephysdir)
%% Ephys analysis for collective behavior experiment
% Analyze d-bat data

svEphys = 'Analysis_TY/Ephys';
mkdir(svEphys)

%% Load data
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

bout_file = dir('Ext_Beh*','bouts*');       load(fullfile(bout_file.folder,bout_file.name));        % flight bout

% bpos_file = dir(fullfile(fileparts(TTlFile.folder),'Ext_Beh*','bpos*'));
% load(fullfile(bpos_file.folder,bpos_file.name));

% Extract session number and date
dirStr = strsplit(fileparts(TTlFile.folder),'\');
DatasetName = dirStr{end-1};        RecDate = dirStr{end};

%% PARAMETERS AND DEFINITIONS
%=== Basic definitions
n_pairs = length(bat_pairs);                                                                                    % Number of bat pairs
bat_ids = [1:n_tags]';                                                                                          % Bat identities
oth_bat = bat_ids(bat_ids~=imp_bat);                                                                            % Other bats
n_cells = length(TT_unit);                                                                                      % Number of cells
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];                                                                           % Override Room boundaries
nsp_clr = [0.2 0.2 0.2; bat_clr(imp_bat,:)];                                                                    % Colormap for Self VS Other
bflying = logical(bflying);                                                                                     % Convert to logical
bat_dist = zeros(T,n_pairs);                                                                                    % Initialize distance between bat pairs
for i = 1:n_pairs,    bat_dist(:,i) = vecnorm(r(:,:,bat_pairs(i,1))-r(:,:,bat_pairs(i,2)),2,2); end             % Distance between bar pairs
angle = fillmissing(angle,'previous',1);                                                                        % Fill missing entries of heading angle with previous value
for i = 1:n_tags; for j = 1:3; custom_map(:,j,i) = linspace(1,bat_clr(i,j))'; end; end                          % Custom graded colormap
f_smp = f_smp;

%=== Parameters and options
n_rep = 1000;                                                                                                   % Number of repetitions for shuffling (Place Cells)
n_rep_PI = 100;                                                                                                 % Number of repetitions for shuffling (Proximity Index)
rstr_int = [-3 3];                                                                                              % Interval for looking at rasters
times_to_shift = t(randi([0.5*Fs T-0.5*Fs],1,n_rep));                                                           % Time intervals for circshift
bin_size_2D = 0.15;                                                                                             % Bin Size for 2D Firing Maps
bin_size_1D = 0.15;                                                                                             % Bin Size for 1D Firing Maps
fake_cell = 0;                                                                                                  % If using the last cell as a benchmark
min_flights_with_spikes = 3;                                                                                    % Minimum number of flights with spikes to compute spatial info
min_time_2D_fly = 0.2;                                                                                          % Minimum time for 2D maps (flight)
min_time_2D_rst = 5;                                                                                            % Minimum time for 2D maps (rest)
d_th_PI = 0.27;                                                                                                 % Threshold distance for the Proximity Index calculation
options.savefigures = 1;                                                                                        % Save Figures
fig_count = 1;                                                                                                  % Id of the first figure
pf = [2 1 0.25];                                                                                                % Time intervals for periflight rasters and evaluations [furthest during close]
sgm = round(1*Fs);                                                                                              % Number of samples for the local segments (1s)
n_sgm = floor(T/sgm);                                                                                           % Number of complete segments that fit into the session
d_th = 10;                                                                                                      % Threshold distance for temporary classification
d = [.6 .9];                                                                                                    % Threshold distance for classifying close vs far {Optimal = [.6 .9]}
n_split = 500;                                                                                                  % Number of times the flights are splitted in two halves
use_fixed_window = 0;                                                                                           % If using fixed time window
sliding_t = [-1:0.1:0.5];                                                                                       % Left edges for the sliding window, social vs non-social flights
window_s = 0.5;                                                                                                 % Sliding window duration
mdg_lim = 0.03;                                                                                                 % Maximal Deviation from g for low-mobility flights
bonferroni_for_window = 1;                                                                                      % If applying Bonferroni correction (a priori) for the number of tested sliding windows (social nonsocial flights)

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

%=== Check if some units are still available
if isempty(TT_unit)
    disp('No good units for this session');
    return
end

%=== Load and align c3d data, not used at the moment
[~,c3d_TTL_samples] = pulsewidth(AnalogSignals(:,2));
c3d_raw = Markers(round(c3d_TTL_samples(1)):round(c3d_TTL_samples(end)),:,:)./1000;
c3d_raw(c3d_raw==0) = nan;
c3d_raw = squeeze(median(c3d_raw,2,'omitnan'));
c3d_raw(1,:) = r(1,:,imp_bat);
c3d_raw(end,:) = r(end,:,imp_bat);
c3d_t = [0:size(c3d_raw,1)-1]'./AnalogFrameRate;

%=== Create analysis folder for storing the results
if ~exist('folder_name')
    analysis_directory=fullfile(pwd,['TT&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
else
    current_directory = pwd;
    analysis_directory = [replace(current_directory,'Ephys_Restrictive_Raw',folder_name),'\TT&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
    if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
end




%% GENERATE THE S CELL WITH SPIKES, SURROGATE SPIKES AND AUTOCORRELATION
s = cell(n_cells,n_rep+1);      % Cell with spike times (s) between first and last TTL of the fly session + surrogates
ac = zeros(n_cells,601);        % Autocorrelation function
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

 %% Flight Analysis
clear AllFlight
for nb = 1:n_tags
    trg_bat = nb;  
    
    AllFlight(nb).bat = trg_bat;
    AllFlight(nb).dataset = DatasetName;
    AllFlight(nb).rec = RecDate;
    

    [btakeoff,blanding] = findchanges(bpos,[-1,1:6],1);
    AllFlight(nb).takeoff = find(btakeoff(:,trg_bat)); AllFlight(nb).landing = find(blanding(:,trg_bat)); % Implanted bat takeoffs/landings
    if length(AllFlight(nb).takeoff) > length(AllFlight(nb).landing)
        AllFlight(nb).takeoff = AllFlight(nb).takeoff(1:end-1); % Remove the final takeoff (which does not have the corresponding landing)
    end
    AllFlight(nb).fllen = AllFlight(nb).landing - AllFlight(nb).takeoff + 1; % Flight length
    AllFlight(nb).flnum = length(AllFlight(nb).takeoff);
    AllFlight(nb).longestfl = max(AllFlight(nb).fllen); % longest flight
    AllFlight(nb).shortestfl = min(AllFlight(nb).fllen); % shortest flight
    
    % Sort with flight length
    [fllen_sorted,idx_sorted] = sort(AllFlight(nb).fllen);
    AllFlight(nb).takeoff = AllFlight(nb).takeoff(idx_sorted);
    AllFlight(nb).landing = AllFlight(nb).landing(idx_sorted);
    AllFlight(nb).fllen = AllFlight(nb).fllen(idx_sorted);
    
    % %===FIGURE: 
    % figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.3 0.55 0.35]);
    % tl = tiledlayout(1, 3, 'TileSpacing', 'tight');
    % nexttile
    % plot(t(2:end),cumsum(histcounts(trgbat_fl.takeoff,1:T)),'k')
    % title('Takeoff of implanted bat')
    % nexttile
    % plot(t(2:end),cumsum(histcounts(trgbat_fl.landing,1:T)),'k')
    % title('Landing of implanted bat')
    % nexttile
    % scatter(trgbat_fl.takeoff/Fs,trgbat_fl.fllen,'k')
    % title('Length of flight')
    
    % Flight interaction and spikes
    pre = 5; % pre-takeoff (s)
    post = 5; % post-landing (s)
    maxfl = AllFlight(nb).longestfl; % maximum length of flight
    len_interp = 500; % length of flight after interporation (timestep)
    alldur = (pre+post)*Fs + len_interp; % Total length of array
    
    othbats = [1:trg_bat-1,trg_bat+1:n_tags];
    
    oth_dist = bat_dist(:,any(bat_pairs == trg_bat,2)); % Distance between implanted bat and others
    oth_pos = bpos(:,othbats); % Position of others 
    oth_fly = double(bflying(:,othbats)); % Flying of others
    oth_r = r(:,:,othbats); % location of other bats
    self_r = r(:,:,trg_bat);
    
    AllFlight(nb).fllim = zeros(AllFlight(nb).flnum,alldur); % 1 when takeoff/landing of imp bat
    AllFlight(nb).othfly = zeros(AllFlight(nb).flnum,alldur,n_tags-1);
    AllFlight(nb).anyothfly = zeros(AllFlight(nb).flnum,alldur);
    AllFlight(nb).numothfly = zeros(AllFlight(nb).flnum,alldur);
    AllFlight(nb).othpos = zeros(AllFlight(nb).flnum,alldur,n_tags-1);
    AllFlight(nb).othdist = zeros(AllFlight(nb).flnum,alldur,n_tags-1);

    AllFlight(nb).othx = zeros(AllFlight(nb).flnum,alldur,n_tags-1);
    AllFlight(nb).othy = zeros(AllFlight(nb).flnum,alldur,n_tags-1);
    AllFlight(nb).othz = zeros(AllFlight(nb).flnum,alldur,n_tags-1);

    AllFlight(nb).x = zeros(AllFlight(nb).flnum,alldur);
    AllFlight(nb).y = zeros(AllFlight(nb).flnum,alldur);
    AllFlight(nb).z = zeros(AllFlight(nb).flnum,alldur);
    
    
    for fl = 1:AllFlight(nb).flnum
        idx_takeoff = AllFlight(nb).takeoff(fl);
        idx_landing = AllFlight(nb).landing(fl);
        fllen_now = AllFlight(nb).fllen(fl);
        ttllen_now = (pre+post)*Fs+fllen_now;
    
        idx_pre = idx_takeoff-pre*Fs:idx_takeoff-1;
        idx_post = idx_landing+1:idx_landing+post*Fs; % peri-flight index

        if min(idx_pre) >= 1 && max(idx_post) <= T
    
            AllFlight(nb).fllim(pre*Fs+1) = 1; % takeoff of impbat
            AllFlight(nb).fllim(pre*Fs+len_interp) = 1; % landing of impbat
        
            AllFlight(nb).othfly(fl,1:pre*Fs,:) = oth_fly(idx_pre,:);
            AllFlight(nb).othfly(fl,pre*Fs+len_interp+1:end,:) = oth_fly(idx_post,:); % flying of other bats during imp bat flights
            AllFlight(nb).othfly(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_fly(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllFlight(nb).anyothfly(fl,1:pre*Fs) = any(oth_fly(idx_pre,:),2); AllFlight(nb).anyothfly(fl,pre*Fs+len_interp+1:end) = any(oth_fly(idx_post,:),2); % flying of other bats during imp bat flights
            AllFlight(nb).anyothfly(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(double(any(oth_fly(idx_takeoff:idx_landing,:),2)),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllFlight(nb).numothfly(fl,1:pre*Fs) = sum(oth_fly(idx_pre,:),2); AllFlight(nb).numothfly(fl,pre*Fs+len_interp+1:end) = sum(oth_fly(idx_post,:),2); % flying of other bats during imp bat flights
            AllFlight(nb).numothfly(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(sum(oth_fly(idx_takeoff:idx_landing,:),2),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllFlight(nb).othpos(fl,1:pre*Fs,:) = oth_pos(idx_pre,:); AllFlight(nb).othpos(fl,pre*Fs+len_interp+1:end,:) = oth_pos(idx_post,:); % flying of other bats during imp bat flights
            AllFlight(nb).othpos(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_pos(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllFlight(nb).othdist(fl,1:pre*Fs,:) = oth_dist(idx_pre,:); AllFlight(nb).othdist(fl,pre*Fs+len_interp+1:end,:) = oth_dist(idx_post,:); % flying of other bats during imp bat flights
            AllFlight(nb).othdist(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_dist(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
                
            AllFlight(nb).othx(fl,1:pre*Fs,:) = oth_r(idx_pre,1,:); AllFlight(nb).othx(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,1,:); % location of other bats during flights
            AllFlight(nb).othx(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,1,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
            AllFlight(nb).othy(fl,1:pre*Fs,:) = oth_r(idx_pre,2,:); AllFlight(nb).othy(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,2,:); % location of other bats during flights
            AllFlight(nb).othy(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,2,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
            AllFlight(nb).othz(fl,1:pre*Fs,:) = oth_r(idx_pre,3,:); AllFlight(nb).othz(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,3,:); % location of other bats during flights
            AllFlight(nb).othz(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,3,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
            
            AllFlight(nb).x(fl,1:pre*Fs)                        = self_r(idx_pre,1);
            AllFlight(nb).x(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,1); % Location during flights
            AllFlight(nb).x(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,1),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
            AllFlight(nb).y(fl,1:pre*Fs)                        = self_r(idx_pre,2);
            AllFlight(nb).y(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,2); % Location during flights
            AllFlight(nb).y(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,2),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
            AllFlight(nb).z(fl,1:pre*Fs)                        = self_r(idx_pre,3);
            AllFlight(nb).z(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,3); % Location during flights
            AllFlight(nb).z(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,3),1:(fllen_now-1)/(len_interp-1):fllen_now);
    
        %     trgbat_fl.anyothfly(fl,1:ttllen_now) = any(oth_fly(idx_pre,:),2); % flying of other bats during imp bat flights
        %     trgbat_fl.numothfly(fl,1:ttllen_now) = sum(oth_fly(idx_pre,:),2); % flying of other bats during imp bat flights
        %     trgbat_fl.othpos(fl,1:ttllen_now,:) = oth_pos(idx_pre,:); % Position of other bats during imp bat flights
        %     trgbat_fl.othdist(fl,1:ttllen_now,:) = oth_dist(idx_pre,:); % Distance from other bats during imp bat flights
        end
    end
    
    for nc = 1:n_cells
        AllFlight(nb).spikes(nc).rate = zeros(AllFlight(nb).flnum,alldur);
        AllFlight(nb).spikes(nc).zrate = zeros(AllFlight(nb).flnum,alldur);
        for fl = 1:AllFlight(nb).flnum
            idx_takeoff = AllFlight(nb).takeoff(fl);
            idx_landing = AllFlight(nb).landing(fl);
            fllen_now = AllFlight(nb).fllen(fl);
            idx_pre = idx_takeoff-pre*Fs:idx_takeoff-1;
            idx_post = idx_landing+1:idx_landing+post*Fs; % peri-flight index

            if min(idx_pre) >= 1 && max(idx_post) <= T
    
                AllFlight(nb).spikes(nc).zrate(fl,1:pre*Fs) = zRate(idx_pre,nc);
                AllFlight(nb).spikes(nc).zrate(fl,pre*Fs+len_interp+1:end) = zRate(idx_post,nc);
                AllFlight(nb).spikes(nc).zrate(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(zRate(idx_takeoff:idx_landing,nc),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
                AllFlight(nb).spikes(nc).rate(fl,1:pre*Fs) = Rate(idx_pre,nc);
                AllFlight(nb).spikes(nc).rate(fl,pre*Fs+len_interp+1:end) = Rate(idx_post,nc);
                AllFlight(nb).spikes(nc).rate(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(Rate(idx_takeoff:idx_landing,nc),1:(fllen_now-1)/(len_interp-1):fllen_now);
            end
        end
    end


    % Spike around flights
    %===FIGURE: 
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    tl = tiledlayout(n_cells, 1, 'TileSpacing', 'tight');
    title(tl,sprintf('%d-th bat: Peri-flight spike activity',nb))
    xlabel(tl,'Time')
    ylabel(tl,'Flights')
   
    for nc = 1:n_cells
        nexttile
        imagesc(AllFlight(nb).spikes(nc).zrate)
        
        if nc ~= n_cells
            set(gca,'xtick',[])
        else
            xticks([1,pre*Fs,Fs*(pre+post),alldur])
            xticklabels({[num2str(-pre),' sec'],'takeoff','landing',[num2str(post), 'sec']}) 
        end
        ylabel(['Cell ',num2str(nc)])
        colorbar
    end

    exportgraphics(tl,fullfile('Analysis_TY',['periflight_activity_',num2str(nb),'.png']),'BackgroundColor','white')


    %===FIGURE: Average firing rate
    f = figure(); set(f, 'units', 'normalized', 'outerposition', [0.3 0.3 0.25 0.3]);
    title(sprintf('%d-th bat: Peri-flight spike activity',nb))

    hold on
    for nc = 1:n_cells
        plot(mean(AllFlight(nb).spikes(nc).zrate,1),'DisplayName',['Cell ',num2str(nc)])
    end
    hold off
    xlabel('Time')
    ylabel('Ffiring rate (Hz)')
    xticks([1,pre*Fs,Fs*(pre+post),alldur])
    xticklabels({[num2str(-pre),' sec'],'takeoff','landing',[num2str(post), 'sec']})
    xlim([0 alldur])
    legend('Location','northeastoutside')
    exportgraphics(f,fullfile('Analysis_TY',['periflight_activity_mean_raw',num2str(nb),'.png']),'BackgroundColor','white')

end

close all
save(fullfile(TTlFile.folder,'AllFlight.mat'),'AllFlight')



%% Distance vs firing rate
% dist_rate_corr = zeros(n_tags-1,n_cells); % correlations between distance with other bats and firing rate for implanted bat
% % close all
% 
% % idx_now = [1:pre*Fs];
% idx_now = [pre*Fs+1:pre*Fs+len_interp];
% % idx_now = [pre*Fs+fllen_interp+1;(pre+post)*Fs+fllen_interp];
% % idx_now = [1:(pre+post)*Fs+fllen_interp];
% 
% for nb=1:n_tags-1
%     for nc = 1:n_cells
%         dist_now = AllFlight(imp_bat).othdist(:,idx_now,nb);    
% %         dist_now = dist_now(:);
%         rate_now = AllFlight(imp_bat).spikes(nc).rate(:,idx_now);
% %         rate_now = rate_now(:);
%         corr_now = corr(dist_now(:),rate_now(:));
%         dist_rate_corr(nb,nc) = corr_now;
%         
%         if abs(corr_now) > 0.25
%             fllen_now = AllFlight(imp_bat).flnum;
% %             c = linspace(1,10,n_xq);
%             c = linspace(1,10,len_interp);
% 
% 
%             figure
%             for nf = 1:fllen_now
%                 hold on
%                 scatter(dist_now(nf,:),rate_now(nf,:),3,c,'filled')
%                 title(sprintf('%d-th bat, %d-th cell, corr=%.3f',nb,nc,corr_now))
%                 xlabel('Distance from implanted bat')
%                 ylabel('Zscored firing rate')
%             end
%         end
%         
%     end
% end
% 
% figure
% imagesc(dist_rate_corr)
% title('Distance-firing rate correlation')
% xlabel('Cell')
% ylabel('Bat')
% yticks(1:n_tags)
% colorbar


%% Distance vs firing rate v2
dist_rate_corr = zeros(n_tags-1,n_cells); % correlations between distance with other bats and firing rate for implanted bat
close all

% idx_now = [1:pre*Fs];
idx_now = [pre*Fs+1:pre*Fs+len_interp];
% idx_now = [pre*Fs+fllen_interp+1;(pre+post)*Fs+fllen_interp];
% idx_now = [1:(pre+post)*Fs+fllen_interp];

n_xq = 100;
w_interp = 1 + (len_interp-1) .* (0.5 + 0.5 .* (linspace(-n_xq/2,n_xq/2,n_xq)).^3 / (n_xq/2)^3);

fig_cnt = 0;
for nb=1:n_tags-1
    for nc = 1:n_cells
        dist_now = interp1(AllFlight(imp_bat).othdist(:,idx_now,nb)',w_interp)';
        rate_now = interp1(AllFlight(imp_bat).spikes(nc).zrate(:,idx_now)',w_interp)';
        corr_now = corr(dist_now(:),rate_now(:));
        dist_rate_corr(nb,nc) = corr_now;
        
        if abs(corr_now) > 0.25
            fllen_now = AllFlight(imp_bat).flnum;
            c = linspace(1,10,n_xq);

            f = figure;
            for nf = 1:fllen_now
                hold on
                scatter(dist_now(nf,:),rate_now(nf,:),3,c,'filled')
                title(sprintf('%d-th bat, %d-th cell, corr=%.3f',nb,nc,corr_now))
                xlabel('Distance from implanted bat')
                ylabel('Zscored firing rate')
            end
            hold off
            exportgraphics(f,fullfile('Analysis_TY',['Dist_rate_corr_',num2str(fig_cnt),'.png']),'BackgroundColor','white')
            fig_cnt = fig_cnt + 1;
        end
        
    end
end

f = figure;
imagesc(dist_rate_corr)
title('Distance-firing rate correlation')
xlabel('Cell')
ylabel('Bat')
yticks(1:n_tags)
colorbar

exportgraphics(f,fullfile('Analysis_TY',['Dist_rate_corr_',num2str(fig_cnt),'.png']),'BackgroundColor','none')
close all

%% Other information for trgbat_fl

% Place activity for location of other bats during flight of implanted bat
fig_cnt = 1;
for nb = 1:n_tags-1
    for nc = 1:n_cells
        f = figure;
        set(f, 'units','normalized','outerposition',[0.2 0.25 0.7 0.35]);
        tl = tiledlayout(1,3,'TileSpacing','tight');
        title(tl,sprintf('Spike activity of %d-th cell vs location of %d-th bat during implanted bat flights',nc,nb))

        x = AllFlight(imp_bat).othx(:,:,nb);
        y = AllFlight(imp_bat).othy(:,:,nb);
        z = AllFlight(imp_bat).othz(:,:,nb);
        c = AllFlight(imp_bat).spikes(nc).rate;
        c = c(:);
        c = c/max(c);
        c = 100.^c;
    
        nexttile
        scatter(x(:),y(:),5,c,'filled','MarkerFaceAlpha',0.15)
        nexttile
        scatter(y(:),z(:),3,c,'filled','MarkerFaceAlpha',0.15)
        nexttile
        scatter(x(:),z(:),3,c,'filled','MarkerFaceAlpha',0.15)
        xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
        colormap hot

        exportgraphics(f,fullfile('Analysis_TY',['PlaceMap_',num2str(fig_cnt),'.png']),'BackgroundColor','white')
        fig_cnt = fig_cnt + 1;
        close gcf
        pause(0.1)
    end
end
%%

% Place activity for location of others/self during flight of others/implanted bats
for nb = 1:n_tags
    for nc = 1:n_cells
        f = figure;
        set(f, 'units','normalized','outerposition',[0.2 0.25 0.7 0.35]);
        tl = tiledlayout(1,3,'TileSpacing','tight');
        title(tl,sprintf('Spike activity of %d-th cell vs location of %d-th bat during %d-th bat flights',nc,nb,nb))
        
        x = AllFlight(nb).x;
        y = AllFlight(nb).y;
        z = AllFlight(nb).z;
        c = AllFlight(nb).spikes(nc).rate;
        c = c(:);
        c = c/max(c);
        c = 100.^c;
        
        nexttile
        scatter(x(:),y(:),5,c,'filled','MarkerFaceAlpha',0.15)
        nexttile
        scatter(y(:),z(:),3,c,'filled','MarkerFaceAlpha',0.15)
        nexttile
        scatter(x(:),z(:),3,c,'filled','MarkerFaceAlpha',0.15)
        xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
        colormap hot

        exportgraphics(f,fullfile('Analysis_TY',['PlaceMap_',num2str(fig_cnt),'.png']),'BackgroundColor','white')
        fig_cnt = fig_cnt + 1;
        close gcf
        pause(0.1)
    end
end


%% Bout Analysis
clear AllBout
for nb = 1:n_tags
    trg_bat = nb;  
    
    AllBout(nb).bat = trg_bat;
    AllBout(nb).dataset = DatasetName;
    AllBout(nb).rec = RecDate;
    

    [btakeoff,blanding] = findchanges(bout,0,1);
    AllBout(nb).boutstart = find(btakeoff(:,trg_bat)); AllBout(nb).boutend = find(blanding(:,trg_bat)); % Implanted bat takeoffs/landings
    if length(AllBout(nb).boutstart) > length(AllBout(nb).boutend)
        AllBout(nb).boutstart = AllBout(nb).boutstart(1:end-1); % Remove the final takeoff (which does not have the corresponding landing)
    end

    AllBout(nb).boutlen = AllBout(nb).boutend - AllBout(nb).boutstart + 1; % Flight length
    AllBout(nb).boutnum = length(AllBout(nb).boutstart);
    AllBout(nb).longestbout = max(AllBout(nb).boutlen); % longest flight
    AllBout(nb).shortestfl = min(AllBout(nb).boutlen); % shortest flight
    
    % Sort with flight length
    [boutlen_sorted,idx_sorted] = sort(AllBout(nb).boutlen);
    AllBout(nb).boutstart = AllBout(nb).boutstart(idx_sorted);
    AllBout(nb).boutend = AllBout(nb).boutend(idx_sorted);
    AllBout(nb).boutlen = AllBout(nb).boutlen(idx_sorted);
    
    % %===FIGURE: 
    % figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.3 0.55 0.35]);
    % tl = tiledlayout(1, 3, 'TileSpacing', 'tight');
    % nexttile
    % plot(t(2:end),cumsum(histcounts(trgbat_fl.takeoff,1:T)),'k')
    % title('Takeoff of implanted bat')
    % nexttile
    % plot(t(2:end),cumsum(histcounts(trgbat_fl.landing,1:T)),'k')
    % title('Landing of implanted bat')
    % nexttile
    % scatter(trgbat_fl.takeoff/Fs,trgbat_fl.fllen,'k')
    % title('Length of flight')
    
    % Flight interaction and spikes
    pre = 10; % pre-takeoff (s)
    post = 10; % post-landing (s)
    len_interp = 1000; % length of flight after interporation (timestep)
    alldur = (pre+post)*Fs + len_interp; % Total length of array
    
    othbats = [1:trg_bat-1,trg_bat+1:n_tags];
    
%     oth_dist = bat_dist(:,any(bat_pairs == trg_bat,2)); % Distance between implanted bat and others
%     oth_pos = bpos(:,othbats); % Position of others 
    oth_bout = double(bout(:,othbats)); % Flying of others
    oth_r = r(:,:,othbats); % location of other bats
    self_r = r(:,:,trg_bat);
    
    AllBout(nb).boutlim = zeros(AllBout(nb).boutnum,alldur); % 1 when takeoff/landing of imp bat
    AllBout(nb).othbout = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
    AllBout(nb).anyothbout = zeros(AllBout(nb).boutnum,alldur);
    AllBout(nb).numothbout = zeros(AllBout(nb).boutnum,alldur);
%     AllBout(nb).othpos = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
%     AllBout(nb).othdist = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
% 
%     AllBout(nb).othx = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
%     AllBout(nb).othy = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
%     AllBout(nb).othz = zeros(AllBout(nb).boutnum,alldur,n_tags-1);
% 
%     AllBout(nb).x = zeros(AllBout(nb).boutnum,alldur);
%     AllBout(nb).y = zeros(AllBout(nb).boutnum,alldur);
%     AllBout(nb).z = zeros(AllBout(nb).boutnum,alldur);
    
    
    for fl = 1:AllBout(nb).boutnum
        idx_takeoff = AllBout(nb).boutstart(fl);
        idx_landing = AllBout(nb).boutend(fl);
        fllen_now = AllBout(nb).boutlen(fl);
        ttllen_now = (pre+post)*Fs+fllen_now;
    
        idx_pre = idx_takeoff-pre*Fs:idx_takeoff-1;
        idx_post = idx_landing+1:idx_landing+post*Fs; % peri-flight index

        if min(idx_pre) >= 1 && max(idx_post) <= T
    
            AllBout(nb).boutlim(pre*Fs+1) = 1; % takeoff of impbat
            AllBout(nb).boutlim(pre*Fs+len_interp) = 1; % landing of impbat
        
            AllBout(nb).othbout(fl,1:pre*Fs,:) = oth_bout(idx_pre,:); AllBout(nb).othbout(fl,pre*Fs+len_interp+1:end,:) = oth_bout(idx_post,:); % flying of other bats during imp bat flights
            AllBout(nb).othbout(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_bout(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllBout(nb).anyothbout(fl,1:pre*Fs) = any(oth_bout(idx_pre,:),2); AllBout(nb).anyothbout(fl,pre*Fs+len_interp+1:end) = any(oth_bout(idx_post,:),2); % flying of other bats during imp bat flights
            AllBout(nb).anyothbout(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(double(any(oth_bout(idx_takeoff:idx_landing,:),2)),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
            AllBout(nb).numothbout(fl,1:pre*Fs) = sum(oth_bout(idx_pre,:),2); AllBout(nb).numothbout(fl,pre*Fs+len_interp+1:end) = sum(oth_bout(idx_post,:),2); % flying of other bats during imp bat flights
            AllBout(nb).numothbout(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(sum(oth_bout(idx_takeoff:idx_landing,:),2),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
%             AllBout(nb).othpos(fl,1:pre*Fs,:) = oth_pos(idx_pre,:); AllBout(nb).othpos(fl,pre*Fs+len_interp+1:end,:) = oth_pos(idx_post,:); % flying of other bats during imp bat flights
%             AllBout(nb).othpos(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_pos(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
%         
%             AllBout(nb).othdist(fl,1:pre*Fs,:) = oth_dist(idx_pre,:); AllBout(nb).othdist(fl,pre*Fs+len_interp+1:end,:) = oth_dist(idx_post,:); % flying of other bats during imp bat flights
%             AllBout(nb).othdist(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_dist(idx_takeoff:idx_landing,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
%                 
%             AllBout(nb).othx(fl,1:pre*Fs,:) = oth_r(idx_pre,1,:); AllBout(nb).othx(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,1,:); % location of other bats during flights
%             AllBout(nb).othx(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,1,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
%     
%             AllBout(nb).othy(fl,1:pre*Fs,:) = oth_r(idx_pre,2,:); AllBout(nb).othy(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,2,:); % location of other bats during flights
%             AllBout(nb).othy(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,2,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
%     
%             AllBout(nb).othz(fl,1:pre*Fs,:) = oth_r(idx_pre,3,:); AllBout(nb).othz(fl,pre*Fs+len_interp+1:end,:) = oth_r(idx_post,3,:); % location of other bats during flights
%             AllBout(nb).othz(fl,pre*Fs+1:pre*Fs+len_interp,:) = interp1(oth_r(idx_takeoff:idx_landing,3,:),1:(fllen_now-1)/(len_interp-1):fllen_now);
%     
%             
%             AllBout(nb).x(fl,1:pre*Fs)                        = self_r(idx_pre,1);
%             AllBout(nb).x(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,1); % Location during flights
%             AllBout(nb).x(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,1),1:(fllen_now-1)/(len_interp-1):fllen_now);
%     
%             AllBout(nb).y(fl,1:pre*Fs)                        = self_r(idx_pre,2);
%             AllBout(nb).y(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,2); % Location during flights
%             AllBout(nb).y(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,2),1:(fllen_now-1)/(len_interp-1):fllen_now);
%     
%             AllBout(nb).z(fl,1:pre*Fs)                        = self_r(idx_pre,3);
%             AllBout(nb).z(fl,pre*Fs+len_interp+1:end)       = self_r(idx_post,3); % Location during flights
%             AllBout(nb).z(fl,pre*Fs+1:pre*Fs+len_interp)    = interp1(self_r(idx_takeoff:idx_landing,3),1:(fllen_now-1)/(len_interp-1):fllen_now);

       
    %     trgbat_fl.anyothfly(fl,1:ttllen_now) = any(oth_fly(idx_pre,:),2); % flying of other bats during imp bat flights
    %     trgbat_fl.numothfly(fl,1:ttllen_now) = sum(oth_fly(idx_pre,:),2); % flying of other bats during imp bat flights
    %     trgbat_fl.othpos(fl,1:ttllen_now,:) = oth_pos(idx_pre,:); % Position of other bats during imp bat flights
    %     trgbat_fl.othdist(fl,1:ttllen_now,:) = oth_dist(idx_pre,:); % Distance from other bats during imp bat flights
        end
    end
    
    for nc = 1:n_cells
        AllBout(nb).spikes(nc).rate = zeros(AllBout(nb).boutnum,alldur);
        AllBout(nb).spikes(nc).zrate = zeros(AllBout(nb).boutnum,alldur);
        for fl = 1:AllBout(nb).boutnum
            idx_takeoff = AllBout(nb).boutstart(fl);
            idx_landing = AllBout(nb).boutend(fl);
            fllen_now = AllBout(nb).boutlen(fl);
            idx_pre = idx_takeoff-pre*Fs:idx_takeoff-1;
            idx_post = idx_landing+1:idx_landing+post*Fs; % peri-flight index
    
            if min(idx_pre) >= 1 && max(idx_post) <= T
                AllBout(nb).spikes(nc).zrate(fl,1:pre*Fs) = zRate(idx_pre,nc);
                AllBout(nb).spikes(nc).zrate(fl,pre*Fs+len_interp+1:end) = zRate(idx_post,nc);
                AllBout(nb).spikes(nc).zrate(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(zRate(idx_takeoff:idx_landing,nc),1:(fllen_now-1)/(len_interp-1):fllen_now);
        
                AllBout(nb).spikes(nc).rate(fl,1:pre*Fs) = Rate(idx_pre,nc);
                AllBout(nb).spikes(nc).rate(fl,pre*Fs+len_interp+1:end) = Rate(idx_post,nc);
                AllBout(nb).spikes(nc).rate(fl,pre*Fs+1:pre*Fs+len_interp) = interp1(Rate(idx_takeoff:idx_landing,nc),1:(fllen_now-1)/(len_interp-1):fllen_now);
            end
        end
    end


    % Spike around flights
    %===FIGURE: 
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    tl = tiledlayout(n_cells, 1, 'TileSpacing', 'tight');
    title(tl,sprintf('%d-th bat: Peri-bout spike activity',nb))
    xlabel(tl,'Time')
    ylabel(tl,'Flights')
   
    for nc = 1:n_cells
        nexttile
        imagesc(AllBout(nb).spikes(nc).zrate)
        
        if nc ~= n_cells
            set(gca,'xtick',[])
        else
            xticks([1,pre*Fs,Fs*(pre+post),alldur])
            xticklabels({[num2str(-pre),' sec'],'Bout start','Bout end',[num2str(post), 'sec']}) 
        end
        ylabel(['Cell ',num2str(nc)])
        colorbar
    end

    exportgraphics(tl,fullfile('Analysis_TY',['peribout_activity_',num2str(nb),'.png']),'BackgroundColor','white')


    %===FIGURE: Average firing rate
    f = figure(); set(f, 'units', 'normalized', 'outerposition', [0.3 0.3 0.25 0.3]);
    title(sprintf('%d-th bat: Peri-bout spike activity',nb))

    hold on
    for nc = 1:n_cells
        plot(mean(AllBout(nb).spikes(nc).zrate,1),'DisplayName',['Cell ',num2str(nc)])
    end
    hold off
    xlabel('Time')
    ylabel('Ffiring rate (Hz)')
    xticks([1,pre*Fs,Fs*(pre+post),alldur])
    xticklabels({[num2str(-pre),' sec'],'takeoff','landing',[num2str(post), 'sec']})
    xlim([0 len_interp+(pre+post)*Fs])
    legend('Location','northeastoutside')
    exportgraphics(f,fullfile('Analysis_TY',['peribout_activity_mean_raw',num2str(nb),'.png']),'BackgroundColor','white')

end

close all
save(fullfile(TTlFile.folder,'AllBout.mat'),'AllBout')

%% Firing on trajectory
% figure


cd(currdir)