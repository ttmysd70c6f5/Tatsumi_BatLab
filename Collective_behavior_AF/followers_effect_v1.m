function followers_effect()
% Analysis of folloers' effects
% Concatenate all sessions to increase the number of flights

%% Directory
currdir = split(cd,'\');
behdir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Behavior',currdir{end}); % behavioral data
ephydir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Ephys\Ephys_Restrictive_Raw',currdir{end}); % Ephys data

sessions = dir(fullfile(behdir,'**/Extracted_Behavior*.mat')); % extract paths of sessions
TTLs = dir(fullfile(ephydir,'**/Ext_Ephys_1/*c3d_1-Bat_Cluster.mat')); % reward signals

load(fullfile(cd, 'Extracted_All.mat')) % common variables across sessions

%% Options for figures
% Parameters for IFI plot
figopt.hedges = 0:1:60; % edges of IFI histogram (sec)
figopt.ifixlim = 60; % sec
figopt.loghedges = 0:0.1:4; % edges of IFI histogram (log10 sec)
figopt.logifixlim = 4; % log10sec

% Parameters for IFI return plot
figopt.ifimax = 60; % sec
figopt.ifibin = 1; % sec
figopt.logifimax = 4; % log10sec
figopt.logifibin = 0.1; % log10sec

% Parameters for CC
figopt.cc_bin = 0.5; % bin duration for cross-correlogram (sec)
figopt.cc_maxlag = 30; % maximum lag for cross-correlogram (sec)
figopt.cc_width = 2; % window size of mask for convlution (sec)
figopt.cc_mask = ones(1,figopt.cc_width/figopt.cc_bin)/(figopt.cc_width/figopt.cc_bin); % mask for convolution smoothing
figopt.cc_ite = 100; % number of iteration for shift predictor

figopt.wcc_wsize = 10; % window size (csec)
figopt.wcc_winc = 1; % window increment (csec)
figopt.wcc_maxlag = 1000; % maximum lag (csec)
figopt.wcc_laginc = 50; % lag increment (csec)


%% Flight segmentation to bouts
%% Compute and segment IFIs with multiple methods
% seg_alg = {'k-means','GWM','GWM combined'}; % clustering algorithm
% num_alg = length(seg_alg); % number of algorithms
% 
% % Arrays and cells to save results
% IFI_all = cell(1,n_tags); % IFI concatenated for all sessions
% IFI_log_all = cell(1,n_tags); % IFI (log10)
% IFI_day = cell(1,n_tags); % Experiment day
% IFI_t_all = cell(1,n_tags); % Time of the flights corresponding to IFIs
% 
% IFI_combined = []; % IFI combined for all bats
% IFI_log_combined = []; % log10IFI combined for all bats
% IFI_day_combined = []; % Day
% IFI_b_combined = []; % Bat
% 
% IFI_cluster_all = cell(num_alg,n_tags); % cluster identification allocated to each IFI value (1=flight bouts)
% 
% bmode_all = cell(num_alg,1); % binary series of flight bouts (bouts = true) for all sessions
% bout_num_all = zeros(num_alg,n_tags); % number of bouts for all sessions
% bout_len_all = cell(num_alg,n_tags); % length of bouts
% bout_sz_all = cell(num_alg,n_tags); % number of flights in bouts
% IBI_all = cell(num_alg,n_tags); % inter bout intervals (sec)
% IBI_log_all = cell(num_alg,n_tags); % inter bout intervals (log10 sec)
% 
% % Computation of IFI
% for ss = 1:length(sessions)
%     load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying')
%     [btakeoff,~] = findchanges(bflying,0,1);
%       
%     for i=1:n_tags
%         % compute IFI
%         ifi = diff(find(btakeoff(:,i))); % IFI (sec)
%         ifi_log = log10(ifi/Fs); % IFI (log10 sec)
%         ifi_day = repmat(ss,[length(ifi),1]); % Experiment day
%         ifi_t = find(btakeoff(:,i)); ifi_t = ifi_t(1:end-1); % timings of takeoffs
%         ifi_b = repmat(i,length(ifi),1);
%         
%         % concatenate
%         IFI_all{1,i} = [IFI_all{1,i}; ifi];
%         IFI_log_all{1,i} = [IFI_log_all{1,i}; ifi_log];
%         IFI_day{1,i} = [IFI_day{1,i}; ifi_day];
%         IFI_t_all{1,i} = [IFI_t_all{1,i}; ifi_t];
%         
%         IFI_combined = [IFI_combined; ifi];
%         IFI_log_combined = [IFI_log_combined; ifi_log];
%         IFI_day_combined = [IFI_day_combined; ifi_day];
%         IFI_b_combined = [IFI_b_combined; ifi_b];
%     end
% end
% 
% % Flight segmentation
% for alg=1:num_alg
%     alg_now = seg_alg{alg};
%     
%     switch alg_now
%     % Flight segmentation applying clustering algorithm to each bat separately
%     case 'k-means'
%         for i=1:n_tags
%             % k-means
%             ifi_all = IFI_log_all{1,i};
% %             ifi_all = IFI_all{1,i};
%             [ifi_idx,centroid] = kmeans(ifi_all,2); % kmeans clustering
%             if centroid(1) > centroid(2) % Allocate ID=1 for smaller cluster
%                 ifi_idx = 3 - ifi_idx; % Allocate cluster 1 to IFI group with smaller values
%             end
% 
%             IFI_cluster_all{1,i} = ifi_idx;
%         end
% 
%     case 'GWM' % Gaussian mixture model (GWM)
%         rng(2);
%         k = 2; % number of GWM clusters
%         option_gwm = statset('MaxIter',1000); % Options for GWM.
%         for i=1:n_tags    
%             % GWM
% %             ifi_all = IFI_all{1,i};
%             ifi_all = IFI_log_all{1,i};
%             gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal', ...
%                 'SharedCovariance',false,'Options',option_gwm); % Fitted GMM
%             ifi_idx = cluster(gmfit,ifi_all); % Cluster index 
% 
%             if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller cluster
%                 ifi_idx = 3 - ifi_idx;
%             end
% 
%             IFI_cluster_all{2,i} = ifi_idx;
%         end
% 
%     % Flight segmentation applying clustering algorithm to all bats at once
%     case 'GWM combined' % Gaussian mixture model (GWM)
%         rng(2);
%         k = 2; % number of GWM clusters
%         option_gwm = statset('MaxIter',1000); % Options for GWM.
% %         ifi_all = IFI_combined;
%         ifi_all = IFI_log_combined;
%         gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal', ...
%             'SharedCovariance',false,'Options',option_gwm); % Fitted GMM
%         ifi_idx = cluster(gmfit,ifi_all); % Cluster index 
% 
%         if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller clusters
%             ifi_idx = 3 - ifi_idx;
%         end
%         
%         for i=1:n_tags
%             IFI_cluster_all{3,i} = ifi_idx(IFI_b_combined==i);
%         end
%     end
% end
% 
% 
% % Create binary timeseries data for flight bouts
% for ss=1:length(sessions)
%     load(fullfile(sessions(ss).folder, sessions(ss).name),'T')
%     
%     % variables to save results for each sessions
%     bmode = cell(num_alg,1);
%     bout_num = zeros(num_alg,n_tags);
%     bout_len = cell(num_alg,n_tags);
%     bout_sz = cell(num_alg,n_tags);
%     IBI = cell(num_alg,n_tags);
%     IBI_log = cell(num_alg,n_tags);
% 
%     IFI = cell(1,n_tags); % Inter flight intervals
%     IFI_log = cell(1,n_tags); % Inter flight intervals (log10)
%     IFI_t = cell(1,n_tags); % time of IFI
%     IFI_cluster = cell(num_alg,n_tags); % cluster number of IFI    
%     
%     % compute binary series of flight bouts with different algorithms
%     for alg=1:num_alg
%         % initialize
%         if ss==1
%             bmode_all{alg,1} = [];
%         end
%         
%         bmode_now = zeros(T,n_tags);
%         
%         for i=1:n_tags
%             ifi = IFI_all{1,i}(IFI_day{1,i}==ss);
%             ifi_idx = IFI_cluster_all{alg,i}(IFI_day{1,i}==ss);
%             ifi_t = IFI_t_all{1,i}(IFI_day{1,i}==ss);
%             IFI{1,i} = ifi; IFI_log{1,i} = ifi_log; IFI_t{1,i} = ifi_t; IFI_cluster{alg,i} = ifi_idx; % save for each session
% 
%             bout_cnt = 1;
%             bout_sz_now = []; % number of flight in a bout
%             for j=1:length(ifi) % skip the first bout
%                 if ifi_t(find(ifi_idx==1,1,'first')) > 5*60*Fs
%                     if ifi_idx(j) == 1
%                         bmode_now(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
%                         bout_cnt = bout_cnt + 1; % count of flight in bout
%                     else
%                         if bout_cnt > 1
%                             bout_sz_now = [bout_sz_now; bout_cnt];
%                             bout_cnt = 1; % initialize
%                         end
%                     end
%                     
%                     
%                 else            
%                     if ifi_idx(1) == 1 % if the initial IFI belongs to the first bout, skip until the first index after the first bout
%                         if j >= find(ifi_idx==2,1,'first') && ifi_idx(j) == 1
%                             bmode_now(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
%                             bout_cnt = bout_cnt + 1; % count of flight in bout
%                         else
%                             if bout_cnt > 1
%                                 bout_sz_now = [bout_sz_now; bout_cnt];
%                                 bout_cnt = 1; % initialize
%                             end
%                         end
%                         
%                         
%                     elseif ifi_idx(1) == 2 % if the initial IFI does not belong to the first bout, skip after the first bout
%                         if find(ifi_idx==1,1,'first') < find(ifi_idx==2,1,'last')
%                             if j >= find(ifi_idx(find(ifi_idx==1,1,'first'):end)==2,1,'first') && ifi_idx(j) == 1
%                                 bmode_now(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
%                                 bout_cnt = bout_cnt + 1; % count of flight in bout
%                             else
%                                 if bout_cnt > 1
%                                     bout_sz_now = [bout_sz_now; bout_cnt];
%                                     bout_cnt = 1; % initialize
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             if bout_cnt > 1
%                 bout_sz_now = [bout_sz_now; bout_cnt];
%             end
%             
%             % size of bouts
%             bout_sz{alg,i} = bout_sz_now;
%             bout_sz_all{alg,i} = [bout_sz_all{alg,i}; bout_sz_now];
%         end
%         
%         bmode{alg,1} = bmode_now;
%         [btakeoff,~] = findchanges(bmode_now,0,1);
%         bout_num(alg,:) = sum(btakeoff);
%         
%         bmode_all{alg,1} = [bmode_all{alg,1}; bmode_now];
%         bout_num_all(alg,:) = bout_num_all(alg,:) + sum(btakeoff);
%         
%         % Inter-bout interval
%         [btakeoff,~] = findchanges(bmode_now,0,1);
%         for i=1:n_tags
%              ibi = diff(find(btakeoff(:,i)))/Fs;
%              ibi_log = log10(ibi);
%              
%              IBI{alg,i} = ibi;
%              IBI_log{alg,i} = ibi_log;
%              IBI_all{alg,i} = [IBI_all{alg,i}; ibi];
%              IBI_log_all{alg,i} = [IBI_log_all{alg,i}; ibi_log];
%         end
%         
%         % length of bouts
%         for i=1:n_tags
%             [btakeoff, blanding] = findchanges(bmode_now,0,1); % timings of initial flight of bouts
%             bout_len_now = (find(blanding(:,i)) - find(btakeoff(:,i)))/Fs; % length of flight bouts
%             bout_len{alg,i} = bout_len_now;
%             bout_len_all{alg,i} = [bout_len_all{alg,i}; bout_len_now]; % append
%         end
% 
%     end
%     
%     %=== SAVE
%     save(fullfile(sessions(ss).folder,'IFI.mat'),'IFI','IFI_log','IFI_t','IFI_cluster','seg_alg','num_alg')
%     save(fullfile(sessions(ss).folder,'bouts.mat'),'bmode','bout_num','bout_len','bout_sz','IBI','IBI_log','seg_alg','num_alg')
%     clear IFI IFI_log IFI_t IFI_cluster bmode bout_num bout_len bout_sz IBI IBI_log
% end
% 
% %=== SAVE
% save(fullfile(cd,'IFI_all.mat'),'IFI_all','IFI_log_all','IFI_day','IFI_t_all','IFI_cluster_all','IFI_combined','IFI_log_combined','seg_alg','num_alg')
% save(fullfile(cd,'bouts_all.mat'),'bmode_all','bout_num_all','bout_len_all','bout_sz_all','IBI_all','IBI_log_all','seg_alg','num_alg')


%% FIGURE to compare segmentation by different methods
% %=== Figure: clustering of IFI
% load(fullfile(cd, 'IFI_all.mat'))
% 
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% tl = tiledlayout(3, n_tags, 'TileSpacing', 'tight');
% title(tl,sprintf('IFI histogram (all sessions, algorithm=%s,%s,%s)',seg_alg{1},seg_alg{2},seg_alg{3}))
% xlabel(tl,'Log10 time (sec)')
% ylabel(tl,'Counts')
% 
% for alg=1:num_alg    
%     for i=1:n_tags
%         nexttile(); hold on
%         for cl = 1:2
%             ifi = IFI_log_all{1,i}; ifi_idx = IFI_cluster_all{alg,i};
%             histogram(ifi(ifi_idx==cl),0:0.1:4)
%         end
%         hold off
%         title(bat_nms(i,:),'Color',bat_clr(i,:))
%         legend('Bout', 'Inter-bout')
%     end
% end
% saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('IFI_segment_all_compare.png')))
% 
% %=== FIGURE: flight/silent mode segmentation
% for alg=1:num_alg
%     for ss=1:length(sessions)
%         load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying','t','f_num')
%         load(fullfile(sessions(ss).folder, 'bouts.mat'),'bmode','bout_num','seg_alg')
%         bmode_now = bmode{alg,1};
% 
%         figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%         tl = tiledlayout(n_tags,1,'TileSpacing','tight');
%         title(tl,sprintf('Flight mode segmentation (algorithm=%s): day%d',seg_alg{alg},ss))
%         xlabel(tl, 'Time(s)');
%         for i = 1:n_tags
%             ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
%             area(t,bflying(:,i)*5,'FaceAlpha',0.8,'LineStyle','none');  hold on;
%             area(t,bmode_now(:,i)*3,'FaceAlpha',0.3,'LineStyle','-','FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:))
%             hold off
%             title(sprintf('%d flights, %d bouts',f_num(i),bout_num(alg,i)));
%             ylabel(bat_nms(i,:),'Color',bat_clr(i,:))
%         end
%         linkaxes(ax,'x');
%         saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('Mode_segment_%s_day%d.png',seg_alg{alg},ss)))
% 
%         pause(0.3)
%         close gcf
%     end
% end
% 
% 
% %=== FIGURE: number of flight bouts
% for alg=1:num_alg
%     load(fullfile(cd, 'bouts_all.mat'))
%     
%     figure();
%     b = bar(bout_num_all(alg,:),'FaceColor','flat');
%     b.CData = bat_clr;
%     title(sprintf('Number of bouts (all sessions, method=%s)',seg_alg{alg}))
%     xlabel('Bat')
%     ylabel('Counts')
%     xticklabels(cellstr(bat_nms))
% 
%     saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('Flightmode_counts_%s.png',seg_alg{alg})))
%     close gcf
% end
% 
% %=== FIGURE: comparison between different algorithms
% % Convert to table data
% load(fullfile(cd,'bouts_all.mat'))
% tbl = cell(0,0); tbl2 = cell(0,0);
% for alg=1:num_alg
%     for i=1:n_tags
%         bout_num = bout_num_all(alg,i); % number of bouts for all sessions
%         tbl_new = [repmat([seg_alg(alg),bat_nms(i,:)],bout_num,1),num2cell(bout_sz_all{alg,i}),num2cell(bout_len_all{alg,i})];
%         tbl = [tbl;tbl_new];
%         
%         tbl2_new = [repmat([seg_alg(alg),bat_nms(i,:)],length(IBI_all{alg,i}),1),num2cell(IBI_all{alg,i}),num2cell(IBI_log_all{alg,i})];
%         tbl2 = [tbl2;tbl2_new];
%     end
% end
% tbl = cell2table(tbl,'VariableNames',["Algorithm","Bat","BoutSize","BoutLength"]);
% tbl2 = cell2table(tbl2,'VariableNames',["Algorithm","Bat","IBI","Log10IBI"]);
% tbl.Bat = categorical(tbl.Bat);
% tbl2.Bat = categorical(tbl2.Bat);
% 
% 
% plt_ori = 'horizontal';
% figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% if strcmp(plt_ori,'horizontal')
%     tl = tiledlayout(1,4,'TileSpacing','tight');
% elseif strcmp(plt_ori,'vertical')
%     tl = tiledlayout(4,1,'TileSpacing','tight');
% else
%     error('Incorrect orientation specified')
% end 
% 
% % Number of bouts
% vals = zeros(n_tags,num_alg);
% for alg=1:num_alg
%     load(fullfile(cd, 'bouts_all.mat'))
%     bmode_all_now = bmode_all{alg};
%     
%     vals(:,alg) = sum(findchanges(bmode_all_now,0,1));
% end
% 
% 
% nexttile()
% if strcmp(plt_ori,'horizontal')
%     barh(vals,'FaceColor','flat');
%     set(gca,'YDir','reverse')
%     xlabel('Counts'); ylabel('Bat')
%     yticklabels(cellstr(bat_nms))
% elseif strcmp(plt_ori,'vertical')
%     bar(vals,'FaceColor','flat');
%     xlabel('Bat'); ylabel('Counts')
%     xticklabels(cellstr(bat_nms))
% else
%     error('Incorrect orientation specified')
% end 
% title(sprintf('Number of bouts (all sessions)'))
% legend(seg_alg,'Location','northeastoutside')
% 
% 
% 
% % Bout Size
% nexttile()
% boxchart(tbl.Bat,tbl.BoutSize,'GroupByColor',tbl.Algorithm,'Orientation',plt_ori)
% title('Number of flight in a bout','FontSize',12)
% if strcmp(plt_ori,'horizontal')
%     set(gca,'xscale','log','YDir','reverse')
%     ylabel('Bat'); xlabel('Log10 counts')
% elseif strcmp(plt_ori,'vertical')
%     set(gca,'yscale','log')
%     xlabel('Bat'); ylabel('Log10 counts')
% else
%     error('Incorrect orientation specified')
% end 
% legend('Location','northeastoutside')
% 
% % Bout Length
% nexttile()
% boxchart(tbl.Bat,tbl.BoutLength,'GroupByColor',tbl.Algorithm,'Orientation',plt_ori)
% title('Length of bouts','FontSize',12)
% if strcmp(plt_ori,'horizontal')
%     set(gca,'xscale','log','YDir','reverse')
%     ylabel('Bat'); xlabel('Bout length (log10 sec)')
% elseif strcmp(plt_ori,'vertical')
%     set(gca,'yscale','log')
%     xlabel('Bat'); ylabel('Bout length (log10 sec)')
% else
%     error('Incorrect orientation specified')
% end 
% legend('Location','northeastoutside')
% 
% % Inter bout intervals
% nexttile()
% boxchart(tbl2.Bat,tbl2.IBI,'GroupByColor',tbl2.Algorithm,'Orientation',plt_ori)
% title('Inter-bout interval','FontSize',12)
% if strcmp(plt_ori,'horizontal')
%     ylabel('Bat'); xlabel('IBI (log10 sec)')
%     set(gca,'xscale','log','YDir','reverse')
% elseif strcmp(plt_ori,'vertical')
%     xlabel('Bat'); ylabel('IBI (log10 sec)')
%     set(gca,'yscale','log')
% else
%     error('Incorrect orientation specified')
% end 
% legend('Location','northeastoutside')
% 
% saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('FlightBout_compare.png')))

%% Compute and segment IFIs with GWM
% seg_alg = {'GWM'}; % clustering algorithm for bout segmentation
% num_alg = length(seg_alg); % number of algorithms

% Arrays and cells to save results
IFI_AllSessions = cell(1,n_tags); % IFI concatenated for all sessions
IFI_log_all = cell(1,n_tags); % IFI (log10)
IFI_day = cell(1,n_tags); % Experiment day
IFI_t_all = cell(1,n_tags); % Time of the flights corresponding to IFIs

IFI_combined = []; % IFI combined for all bats
IFI_log_combined = []; % log10IFI combined for all bats
IFI_day_combined = []; % Day
IFI_b_combined = []; % Bat

IFI_cluster_all = cell(1,n_tags); % cluster identification allocated to each IFI value (1=flight bouts)

bmode_all = zeros(0,n_tags); % binary series of flight bouts (bouts = true) for all sessions
bout_session = zeros(0,1); % session

bout_num_all = zeros(1,n_tags); % number of bouts for all sessions
bout_len_all = cell(1,n_tags); % length of bouts
bout_sz_all = cell(1,n_tags); % number of flights in bouts
IBI_all = cell(1,n_tags); % inter bout intervals (sec)
IBI_log_all = cell(1,n_tags); % inter bout intervals (log10 sec)


% Segmentation methods
seg_alg = 'GWM';

% Computation of IFI
for ss = 1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying')
    [btakeoff,~] = findchanges(bflying,0,1);
      
    for i=1:n_tags
        % compute IFI
        ifi = diff(find(btakeoff(:,i))); % IFI (sec)
        ifi_log = log10(ifi/Fs); % IFI (log10 sec)
        ifi_day = repmat(ss,[length(ifi),1]); % Experiment day
        ifi_t = find(btakeoff(:,i)); ifi_t = ifi_t(1:end-1); % timings of takeoffs
        ifi_b = repmat(i,length(ifi),1);
        
        % concatenate
        IFI_AllSessions{1,i} = [IFI_AllSessions{1,i}; ifi];
        IFI_log_all{1,i} = [IFI_log_all{1,i}; ifi_log];
        IFI_day{1,i} = [IFI_day{1,i}; ifi_day];
        IFI_t_all{1,i} = [IFI_t_all{1,i}; ifi_t];
        
        IFI_combined = [IFI_combined; ifi];
        IFI_log_combined = [IFI_log_combined; ifi_log];
        IFI_day_combined = [IFI_day_combined; ifi_day];
        IFI_b_combined = [IFI_b_combined; ifi_b];
    end
end

% Flight segmentation
   
switch seg_alg
% Flight segmentation applying clustering algorithm to each bat separately
case 'k-means'
    for i=1:n_tags
        % k-means
        ifi_all = IFI_log_all{1,i};
%             ifi_all = IFI_all{1,i};
        [ifi_idx,centroid] = kmeans(ifi_all,2); % kmeans clustering
        if centroid(1) > centroid(2) % Allocate ID=1 for smaller cluster
            ifi_idx = 3 - ifi_idx; % Allocate cluster 1 to IFI group with smaller values
        end

        IFI_cluster_all{1,i} = ifi_idx;
    end

case 'GWM' % Gaussian mixture model (GWM)
    rng(2);
    k = 2; % number of GWM clusters
    option_gwm = statset('MaxIter',1000); % Options for GWM.
    for i=1:n_tags    
        % GWM
%             ifi_all = IFI_all{1,i};
        ifi_all = IFI_log_all{1,i};
        gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal', ...
            'SharedCovariance',false,'Options',option_gwm); % Fitted GMM
        ifi_idx = cluster(gmfit,ifi_all); % Cluster index 

        if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller cluster
            ifi_idx = 3 - ifi_idx;
        end

        IFI_cluster_all{1,i} = ifi_idx;
    end

% Flight segmentation applying clustering algorithm to all bats at once
case 'GWM combined' % Gaussian mixture model (GWM)
    rng(2);
    k = 2; % number of GWM clusters
    option_gwm = statset('MaxIter',1000); % Options for GWM.
%         ifi_all = IFI_combined;
    ifi_all = IFI_log_combined;
    gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal', ...
        'SharedCovariance',false,'Options',option_gwm); % Fitted GMM
    ifi_idx = cluster(gmfit,ifi_all); % Cluster index 

    if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller clusters
        ifi_idx = 3 - ifi_idx;
    end
    
    for i=1:n_tags
        IFI_cluster_all{1,i} = ifi_idx(IFI_b_combined==i);
    end
end

% Create binary timeseries data for flight bouts
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'T')

    % variables to save results for each sessions
    bmode = zeros(T,n_tags);
    bout_len = cell(1,n_tags);
    bout_sz = cell(1,n_tags);
    IBI = cell(1,n_tags);
    IBI_log = cell(1,n_tags);

    IFI = cell(1,n_tags); % Inter flight intervals
    IFI_log = cell(1,n_tags); % Inter flight intervals (log10)
    IFI_t = cell(1,n_tags); % time of IFI
    IFI_cluster = cell(1,n_tags); % cluster number of IFI    
    
    % compute binary series of flight bouts with different algorithms

    for i=1:n_tags
        ifi = IFI_AllSessions{1,i}(IFI_day{1,i}==ss);
        ifi_idx_now = IFI_cluster_all{1,i}(IFI_day{1,i}==ss);
        ifi_t = IFI_t_all{1,i}(IFI_day{1,i}==ss);
        IFI{1,i} = ifi; IFI_log{1,i} = ifi_log; IFI_t{1,i} = ifi_t; IFI_cluster{1,i} = ifi_idx_now; % save for each session

        bout_cnt = 1;
        bout_sz_now = []; % number of flight in a bout
        for j=1:length(ifi) % skip the first bout
            if ifi_t(find(ifi_idx_now==1,1,'first')) > 5*60*Fs
                if ifi_idx_now(j) == 1
                    bmode(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                    bout_cnt = bout_cnt + 1; % count of flight in bout
                else
                    if bout_cnt > 1
                        bout_sz_now = [bout_sz_now; bout_cnt];
                        bout_cnt = 1; % initialize
                    end
                end
                
                
            else            
                if ifi_idx_now(1) == 1 % if the initial IFI belongs to the first bout, skip until the first index after the first bout
                    if j >= find(ifi_idx_now==2,1,'first') && ifi_idx_now(j) == 1
                        bmode(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                        bout_cnt = bout_cnt + 1; % count of flight in bout
                    else
                        if bout_cnt > 1
                            bout_sz_now = [bout_sz_now; bout_cnt];
                            bout_cnt = 1; % initialize
                        end
                    end
                    
                    
                elseif ifi_idx_now(1) == 2 % if the initial IFI does not belong to the first bout, skip after the first bout
                    if find(ifi_idx_now==1,1,'first') < find(ifi_idx_now==2,1,'last')
                        if j >= find(ifi_idx_now(find(ifi_idx_now==1,1,'first'):end)==2,1,'first') && ifi_idx_now(j) == 1
                            bmode(ifi_t(j):ifi_t(j)+ifi(j)-1,i) = 1; % 1 when bats fly with short intervals (flight mode)
                            bout_cnt = bout_cnt + 1; % count of flight in bout
                        else
                            if bout_cnt > 1
                                bout_sz_now = [bout_sz_now; bout_cnt];
                                bout_cnt = 1; % initialize
                            end
                        end
                    end
                end
            end
        end
        
        if bout_cnt > 1
            bout_sz_now = [bout_sz_now; bout_cnt];
        end
        
        % size of bouts
        bout_sz{1,i} = bout_sz_now;
        bout_sz_all{1,i} = [bout_sz_all{1,i}; bout_sz_now];
    end
    
    [btakeoff,~] = findchanges(bmode,0,1);
    bout_num = sum(btakeoff);
    
    bmode_all = [bmode_all; bmode];
    bout_num_all = bout_num_all + bout_num;
    bout_session = [bout_session; repmat(ss,[T,1])];
    
    % Inter-bout interval
    [btakeoff,~] = findchanges(bmode,0,1);
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
        [btakeoff, blanding] = findchanges(bmode,0,1); % timings of initial flight of bouts
        bout_len_now = (find(blanding(:,i)) - find(btakeoff(:,i)))/Fs; % length of flight bouts
        bout_len{1,i} = bout_len_now;
        bout_len_all{1,i} = [bout_len_all{1,i}; bout_len_now]; % append
    end

    
    %=== SAVE
    save(fullfile(sessions(ss).folder,'IFI.mat'),'IFI','IFI_log','IFI_t','IFI_cluster','seg_alg')
    save(fullfile(sessions(ss).folder,'bouts.mat'),'bmode','bout_num','bout_len','bout_sz','IBI','IBI_log','seg_alg')
    clear IFI IFI_log IFI_t IFI_cluster bmode bout_num bout_len bout_sz IBI IBI_log
end

%=== SAVE
save(fullfile(cd,'IFI_all.mat'),'IFI_AllSessions','IFI_log_all','IFI_day','IFI_t_all','IFI_cluster_all','IFI_combined','IFI_log_combined','seg_alg')
save(fullfile(cd,'bouts_all.mat'),'bmode_all','bout_num_all','bout_len_all','bout_sz_all','bout_session','IBI_all','IBI_log_all','seg_alg')

%% FIGURE to compare segmentation by different methods
%=== Figure: clustering of IFI
load(fullfile(cd, 'IFI_all.mat'))

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.3 0.7 0.35]);
tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
title(tl,sprintf('IFI histogram (all sessions, algorithm=%s)',seg_alg))
xlabel(tl,'Log10 time (sec)')
ylabel(tl,'Counts')

for i=1:n_tags
    nexttile(); hold on
    for cl = 1:2
        ifi = IFI_log_all{1,i}; ifi_idx = IFI_cluster_all{1,i};
        histogram(ifi(ifi_idx==cl),0:0.1:4)
    end
    hold off
    title(bat_nms(i,:),'Color',bat_clr(i,:))
    legend('Bout', 'Inter-bout')
end
saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('IFI_segment_all_%s.png',seg_alg)))

%=== FIGURE: flight/silent mode segmentation
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying','t','f_num')
    load(fullfile(sessions(ss).folder, 'bouts.mat'),'bmode','bout_num','seg_alg')

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,1,'TileSpacing','tight');
    title(tl,sprintf('Flight mode segmentation (algorithm=%s): day%d',seg_alg,ss))
    xlabel(tl, 'Time(s)');
    for i = 1:n_tags
        ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
        area(t,bflying(:,i)*5,'FaceAlpha',0.8,'LineStyle','none');  hold on;
        area(t,bmode(:,i)*3,'FaceAlpha',0.3,'LineStyle','-','FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:))
        hold off
        title(sprintf('%d flights, %d bouts',f_num(i),bout_num(1,i)));
        ylabel(bat_nms(i,:),'Color',bat_clr(i,:))
    end
    linkaxes(ax,'x');
    saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('Mode_segment_%s_day%d.png',seg_alg,ss)))

    pause(0.3)
    close gcf
end

%=== FIGURE: number of flight bouts

load(fullfile(cd, 'bouts_all.mat'))

figure();
b = bar(bout_num_all(1,:),'FaceColor','flat');
b.CData = bat_clr;
title(sprintf('Number of bouts (all sessions, method=%s)',seg_alg))
xlabel('Bat')
ylabel('Counts')
xticklabels(cellstr(bat_nms))

saveas(gcf,fullfile(cd,'Figs','Segmentation',sprintf('Flightmode_counts_%s.png',seg_alg)))
close gcf


%% FIGURE: comparison across sessions
%=== FIGURE: Number of bouts for all sessions
load(fullfile(cd, 'bouts_all.mat'),'bout_num_all')
vals = bout_num_all;

figure(); set(gcf, 'units','normalized','outerposition',[0.35 0.3 0.3 0.45]);
b=bar(vals,'FaceColor','flat');
b.CData = bat_clr;
xlabel('Bat'); ylabel('Counts')
xticklabels(cellstr(bat_nms))
title(sprintf('Number of bouts (all sessions), method=%s',seg_alg))


%=== FIGURE: Number of bouts across sessions
load(fullfile(cd,'Extracted_All.mat'),'bat_nms','n_tags')
tbl = cell(0,0);
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder,'bouts.mat'),'bout_num')
    tbl_new = [cellstr(bat_nms),num2cell(repmat(ss,[n_tags,1])),num2cell(bout_num')];
    tbl = [tbl; tbl_new];
end
tbl = cell2table(tbl,'VariableNames',["Bat","Session","NumBout"]);

figure(); set(gcf, 'units','normalized','outerposition',[0.35 0.3 0.3 0.45]);
b=boxchart(categorical(tbl.Bat),tbl.NumBout);
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
xlabel('Bat'); ylabel('Number of bouts')
title(sprintf('Number of bouts (all sessions), method=%s',seg_alg))


% bout size & bout length
% Convert to table data
load(fullfile(cd,'bouts_all.mat'))

% bout size
tbl = cell(0,0);
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder,'bouts.mat'),'bout_num','bout_len','bout_sz','IBI','IBI_log')
    for i=1:n_tags
        tbl_new = [repmat({bat_nms(i,:),ss},bout_num(i),1),...
            num2cell(bout_sz{1,i}),...
            num2cell(bout_len{1,i}),...
            num2cell([NaN;IBI{1,i}]),...
            num2cell([NaN; IBI_log{1,i}])...
            ];
        tbl = [tbl;tbl_new];
    end
end
tbl = cell2table(tbl,'VariableNames',["Bat","sessions","BoutSize","BoutLength","IBI","logIBI"]);

%=== FIGURE: bout size
figure(); set(gcf, 'units','normalized','outerposition',[0.35 0.3 0.3 0.45]);
b=boxchart(categorical(tbl.Bat),tbl.BoutSize);
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
set(gca,'YScale','log')
xlabel('Bat'); ylabel('Number of flights in a bout (log10)')
title(sprintf('Bout size (all sessions), method=%s',seg_alg))

%=== FIGURE: bout length
figure(); set(gcf, 'units','normalized','outerposition',[0.35 0.3 0.3 0.45]);
b=boxchart(categorical(tbl.Bat),tbl.BoutLength);
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
set(gca,'YScale','log')
xlabel('Bat'); ylabel('Length of bouts (log10 sec)')
title(sprintf('Bout length (all sessions), method=%s',seg_alg))

%=== FIGURE: inter-bout interval (IBI)
figure(); set(gcf, 'units','normalized','outerposition',[0.35 0.3 0.3 0.45]);
b=boxchart(categorical(tbl.Bat),tbl.IBI);
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
set(gca,'YScale','log')
xlabel('Bat'); ylabel('Inter-bout intervals (log10 sec)')
title(sprintf('IBI (all sessions), method=%s',seg_alg))


%% Flight segmentation with GWM

% % Compare covariance structures for flight segmentation using GMM
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% tl = tiledlayout(4, n_tags, 'TileSpacing', 'tight');
% 
% rng(2);
% k = 2; % Number of GMM components
% option_gwm = statset('MaxIter',1000);
% 
% th_gmn = sqrt(chi2inv(0.99,2));
% cnt = 1; % count
% Sigma = {'diagonal','full'}; % Options for covariance matrix type
% nSigma = numel(Sigma);
% 
% SharedCovariance = {true,false}; % Indicator for identical or nonidentical covariance matrices
% SCtext = {'true','false'};
% nSC = numel(SharedCovariance);
% 
% for i = 1:nSigma
%     for j = 1:nSC
%         for bb=1:n_tags
%             ifi_all = IFI_log_all{1,bb};
%             
%             gmfit = fitgmdist(ifi_all,k,'CovarianceType',Sigma{i}, ...
%                 'SharedCovariance',SharedCovariance{j},'Options',option_gwm); % Fitted GMM
%             ifi_idx = cluster(gmfit,ifi_all); % Cluster index 
%             
%             if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller cluster
%                 ifi_idx = 3 - ifi_idx;
%             end
%             
%             nexttile(); hold on
%             histogram(ifi_all(ifi_idx==1),0:0.1:4)
%             histogram(ifi_all(ifi_idx==2),0:0.1:4)
%             hold off
%             title(sprintf('Sigma is %s\nSharedCovariance = %s',Sigma{i},SCtext{j}),'FontSize',8)
%             
%             cnt = cnt + 1;
%         end
%     end
% end
% 
% % flight segmentation with GWM
% rng(2);
% k = 2; % Number of GMM components
% option_gwm = statset('MaxIter',1000);
% th_gmn = sqrt(chi2inv(0.99,2));
% 
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.3 0.9 0.35]);
% tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
% title(tl,sprintf('IFI histogram (all sessions)', ss))
% xlabel(tl,'Log10 time (sec)')
% ylabel(tl,'Log10 counts')
% 
% for bb=1:n_tags
%     ifi_all = IFI_log_all{1,bb};
% 
%     gmfit = fitgmdist(ifi_all,k,'CovarianceType','diagonal', ...
%         'SharedCovariance',false,'Options',option_gwm); % Fitted GMM
%     ifi_idx = cluster(gmfit,ifi_all); % Cluster index 
% 
%     if gmfit.mu(1) > gmfit.mu(2) % Allocate ID=1 for smaller cluster
%         ifi_idx = 3 - ifi_idx;
%     end
% 
%     nexttile(); hold on
%     histogram(ifi_all(ifi_idx==1),0:0.1:4)
%     histogram(ifi_all(ifi_idx==2),0:0.1:4)
%     hold off
% 
%     title(bat_nms(i,:),'Color',bat_clr(i,:))
%     legend('Flight', 'Silent')
% end

%% Inter bouts interval (each bat)
% [btakeoff,~] = findchanges(bmode_all_now,0,1);
% IBI = cell(1,n_tags); % inter-bout interval
% for i=1:n_tags
%     ibi = diff(find(btakeoff(:,i))) / Fs; % intervals (sec)
%     IBI{1,i} = log10(ibi); % logarithmize by 10
% end
% 
% % Parameters for IBI plot
% figopt.logibihedges = 0:0.1:5; % edges of IFI histogram (log10 sec)
% figopt.logibixlim = 5; % log10sec
% 
% %=== FIGURE: inter bout interval
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.55 0.9 0.35]);
% tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
% title(tl, 'Inter-bout interval: all sessions');
% xlabel(tl,'Log10 time (sec)'); ylabel(tl,'Count')
% hedge = figopt.logibihedges; ifi_xlim = figopt.logibixlim;
% 
% for i=1:n_tags
%     nexttile(); histogram(IBI{1,i},hedge);
%     set(gca,'yscale','log')
%     title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
%     xlim([0,ifi_xlim])
%     hold on;
%     plot(repmat(median(IBI{1,i}),[1,2]),ylim,'k--');
%     hold off
%     legend('',sprintf('Median: %.0fsec',10^median(IBI{1,i})));
% end

%% Inter bout interval (all bats)
% [btakeoff,~] = findchanges(bmode_all_now,0,1);
% IBI = log10(diff(find(any(btakeoff,2))) / Fs);
% 
% % Parameters for IBI plot
% figopt.logibihedges = 0:0.1:5; % edges of IFI histogram (log10 sec)
% figopt.logibixlim = 5; % log10sec
% 
% hedge = figopt.logibihedges; ifi_xlim = figopt.logibixlim;
% figure();
% histogram(IBI,hedge)
% set(gca,'yscale','log')
% xlim([0,ifi_xlim])
% hold on;
% plot(repmat(median(IBI),[1,2]),ylim,'k--');
% hold off
% legend('',sprintf('Median: %.0fsec',10^median(IBI)));


%% Computation of PSTH for takeoff/landing flights of bouts/all flights
% Load data
load(fullfile(cd,'bouts_all.mat'))
load(fullfile(cd,'Extracted_All.mat'),'bflying_all')
seg_alg = 'GWM';

% Choose data to compute PSTH
% bmode_all_now = bflying_all; psth_type = 'bflying_takeoff'; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% bmode_all_now = bflying_all; psth_type = 'bflying_landing'; [~,btakeoff] = findchanges(bmode_all_now,0,1);
bmode_all_now = bmode_all; psth_type = 'bout_takeoff'; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% bmode_all_now = bmode_all; psth_type = 'bout_landing'; [~,btakeoff] = findchanges(bmode_all_now,0,1);


% PSTH properties
pre = 20; post = 20; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
% ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
% bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
smoothing = 10; % window size of smoothing (samples)

% Create table for psth
tbl = table();

for i = 1:n_tags
    b_trig = bat_nms(i,:); % trigger bat
    trig = find(btakeoff(:,i)); % trigger timings (samples)

    for tr = 1:length(trig)
        if mod(tr,1000)==0
            fprintf('%s: %d/%d \n',b_trig,tr,length(trig))
        end
        for j=1:n_tags
            if i~=j
                b_resp = bat_nms(j,:); % response bat which PSTH is computed
    
                ts = zeros((pre+post)*Fs+1,1); % trig x time x bat(n-1)
                
                if trig(tr)-pre*Fs > 0 && trig(tr)+post*Fs < size(btakeoff,1)
                    ts(:) =  btakeoff(trig(tr)-pre*Fs:trig(tr)+post*Fs,j); % flight of the response bat arround a trigger
                    flight_times = find(ts); % timing of flights around a trigger (samples)
                    tbl_new = cell2table({b_trig,b_resp,t_all(tr),exp_d(tr),flight_times});
                    tbl = [tbl;tbl_new];
                end
            end
        end
    end
end
tbl.Properties.VariableNames = ["TrigBat","RespBat","t","session","FlightTimes"];
save(fullfile(cd,sprintf('PSTH_table_%s.mat',psth_type)),'tbl','pre','post','smoothing','bin_dur')


%=== FIGURE: PSTH of flight bouts for a bat vs any other bats
figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (data=%s, method=%s)',psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')

for i=1:n_tags
    idx = strcmp(tbl.RespBat,bat_nms(i,:)); % table idx to extract target data
    binsize = bin_dur * 1000; % (msec)
    ntrials = sum(idx); % number of trials
    triallen = (pre+post)*Fs+1; % length of a trial (samples)
    

    flight_times_concatenated = cell2mat(tbl.FlightTimes(idx));
    
    nexttile()
    ax = gca;
    ph = psth(flight_times_concatenated, binsize, Fs, ntrials, triallen, smoothing, ax);
    ph.FaceColor = bat_clr(i,:);
    ph.EdgeColor = 'none';
    xticklabels(xticks/1000-pre)
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_summed_%s.png',psth_type)))

%=== FIGURE: PSTH of flight bouts for every bat pairs
figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (data=%s, method=%s)',psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')

for i=1:n_tags
    for j=1:n_tags
        nexttile()
        if i~=j
            idx = strcmp(tbl.TrigBat,bat_nms(i,:)) & strcmp(tbl.RespBat,bat_nms(j,:)); % table idx to extract target data
            binsize = bin_dur * 1000; % (msec)
            ntrials = sum(idx); % number of trials
            triallen = (pre+post)*Fs+1; % length of a trial (samples)
            
        
            flight_times_concatenated = cell2mat(tbl.FlightTimes(idx));
            
            ax = gca;
            ph = psth(flight_times_concatenated, binsize, Fs, ntrials, triallen, smoothing, ax);
            ph.FaceColor = bat_clr(j,:);
            ph.EdgeColor = 'none';
            xticklabels(xticks/1000-pre)
        end
    end
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_allPairs_%s.png',psth_type)))

%% Shuffling test for PSTH for takeoff/landing flights of bouts/all flights
% Load data
load(fullfile(cd,'bouts_all.mat'))
load(fullfile(cd,'Extracted_All.mat'),'bflying_all')
seg_alg = 'GWM';

% Choose data to compute PSTH
% bmode_all_now = bflying_all; psth_type = 'bflying_takeoff'; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% bmode_all_now = bflying_all; psth_type = 'bflying_landing'; [~,btakeoff] = findchanges(bmode_all_now,0,1);
bmode_all_now = bmode_all; psth_type = 'bout_takeoff'; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% bmode_all_now = bmode_all; psth_type = 'bout_landing'; [~,btakeoff] = findchanges(bmode_all_now,0,1);


% PSTH properties
pre = 20; post = 20; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
% ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
% bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
smoothing = 10; % window size of smoothing (samples)
n_ite = 100; % number of iteration


% Create table for psth
tbl_shuffle = table();

for i = 1:n_tags
    % determine the degree of shift randomly
    shft_dist = [2*prctile(IBI_all{1,i},95) size(btakeoff,1)-2*prctile(IBI_all{1,i},95)]; % minimum and maximum shift
    rng(2)
    shft = randi(int32(shft_dist),[1 n_ite]);

    b_trig = bat_nms(i,:); % trigger bat

    for tr = 1:sum(btakeoff(:,i))
        fprintf('%s: completed %d/%d triggers \n',b_trig,tr,sum(btakeoff(:,i)))

        for j=1:n_tags
            if i~=j
                b_resp = bat_nms(j,:); % response bat which PSTH is computed
              
                flight_times = [];
            
                for ite = 1:n_ite
                    trig = find(circshift(btakeoff(:,i),shft(ite))); % trigger timings (samples)        
                    ts = zeros((pre+post)*Fs+1,1); % trig x time x bat(n-1)
                    
                    if trig(tr)-pre*Fs > 0 && trig(tr)+post*Fs < size(btakeoff,1)
                        ts(:) =  btakeoff(trig(tr)-pre*Fs:trig(tr)+post*Fs,j); % flight of the response bat arround a trigger
                        flight_times = [flight_times; find(ts)]; % timing of flights around a trigger (samples)
                    end
                end
                tbl_new = cell2table({b_trig,b_resp,exp_d(tr),flight_times});
                tbl_shuffle = [tbl_shuffle;tbl_new];
            end
        end
    end
end
tbl_shuffle.Properties.VariableNames = ["TrigBat","RespBat","session","FlightTimes"];
save(fullfile(cd,sprintf('PSTH_table_%s_shuffle.mat',psth_type)),'tbl_shuffle','pre','post','smoothing','bin_dur')
fprintf('%s: Shuffled psth was saved!',psth_type)

%%
%=== FIGURE: PSTH of flight bouts for a bat vs any other bats
figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('Shuffled PSTH of bout (data=%s, method=%s)',psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')

for i=1:n_tags
    idx = strcmp(tbl_shuffle.RespBat,bat_nms(i,:)); % table idx to extract target data
    binsize = bin_dur * 1000; % (msec)
    ntrials = sum(idx)*n_ite; % number of trials
    triallen = (pre+post)*Fs+1; % length of a trial (samples)
    

    flight_times_concatenated = cell2mat(tbl_shuffle.FlightTimes(idx));
    
    nexttile()
    ax = gca;
    ph = psth(flight_times_concatenated, binsize, Fs, ntrials, triallen, smoothing, ax);
    ph.FaceColor = bat_clr(i,:);
    ph.EdgeColor = 'none';
    xticklabels(xticks/1000-pre)
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_summed_shuffled_%s.png',psth_type)))

%=== FIGURE: PSTH of flight bouts for every bat pairs
figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,sprintf('Shuffled PSTH of bouts (data=%s, method=%s)',psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')

for i=1:n_tags
    for j=1:n_tags
        nexttile()
        if i~=j
            idx = strcmp(tbl_shuffle.TrigBat,bat_nms(i,:)) & strcmp(tbl_shuffle.RespBat,bat_nms(j,:)); % table idx to extract target data
            binsize = bin_dur * 1000; % (msec)
            ntrials = sum(idx)*n_ite; % number of trials
            triallen = (pre+post)*Fs+1; % length of a trial (samples)
            
        
            flight_times_concatenated = cell2mat(tbl_shuffle.FlightTimes(idx));
            
            ax = gca;
            ph = psth(flight_times_concatenated, binsize, Fs, ntrials, triallen, smoothing, ax);
            ph.FaceColor = bat_clr(j,:);
            ph.EdgeColor = 'none';
            xticklabels(xticks/1000-pre)
        end
    end
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_allPairs_shuffled_%s.png',psth_type)))


%% Statistical test
load(fullfile(cd,'PSTH_table_bout_takeoff.mat'))
load(fullfile(cd,'PSTH_table_bout_takeoff_shuffle.mat'))

th = 0.05;
pvals = zeros(2,n_tags);

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
for i=1:n_tags
    idx = strcmp(tbl.RespBat,bat_nms(i,:)); % table idx to extract target data
    idx_shuffle = strcmp(tbl_shuffle.RespBat,bat_nms(i,:)); % table idx to extract target data
    binsize = bin_dur * 1000; % (msec)
    ntrials = sum(idx); % number of trials
    ntrials_shuffle = sum(idx_shuffle); % number of trials
    triallen = (pre+post)*Fs+1; % length of a trial (samples)
    

    flight_times_concatenated = cell2mat(tbl.FlightTimes(idx));
    flight_times_shuffle_concatenated = cell2mat(tbl_shuffle.FlightTimes(idx_shuffle));

    %Compute raw psth      
    lastBin = binsize * ceil((triallen-1)*(1000/(Fs*binsize))); % the last bin (msec)
    edges = 0 : binsize : lastBin; % edge of psth (msec)
    x = (mod(flight_times_concatenated-1,triallen)+1)*(1000/Fs); % x axis of psth (msec). The first sample is set to 0-th bin
    psthvals = movmean((histc(x,edges)*1000) / (ntrials*binsize), smoothing); % Calculate firing rate (Fz). 1000 is multiplied to convert (1/msec) to (1/sec)

    %Compute shift predictor
    lastBin = binsize * ceil((triallen-1)*(1000/(Fs*binsize))); % the last bin (msec)
    edges = 0 : binsize : lastBin; % edge of psth (msec)
    x = (mod(flight_times_shuffle_concatenated-1,triallen)+1)*(1000/Fs); % x axis of psth (msec). The first sample is set to 0-th bin
    psthvals_shuffle = movmean((histc(x,edges)*1000) / (n_ite*ntrials_shuffle*binsize), smoothing); % Calculate firing rate (Fz). 1000 is multiplied to convert (1/msec) to (1/sec)


    
    [pks,pklocs] = findpeaks(psthvals,'NPeaks',1,'SortStr','descend');
    [troughs,trlocs] = findpeaks(-psthvals,'NPeaks',1,'SortStr','descend');
    lambda = psthvals_shuffle(pklocs);
    pvals(1,i) = min(poisscdf(pks,lambda),1-poisscdf(pks,lambda));
    lambda = psthvals_shuffle(trlocs);
    pvals(2,i) = min(poisscdf(-troughs,lambda),1-poisscdf(-troughs,lambda));

    % FIGURE  
    nexttile()
    ax = gca;
%     bar(psthvals,'EdgeColor','none','FaceColor',bat_clr(i,:));
    h(1) = plot(psthvals,'k','DisplayName','Row');
    hold on
    h(2) = plot(psthvals_shuffle,'k--','DisplayName','Shifted');
    h(3) = plot([pklocs pklocs],ylim,'r--','DisplayName','');
    h(4) = plot([trlocs trlocs],ylim,'b--','DisplayName','');
    hold off
    legend(h(1:2))
    xlim([0 (pre+post)/bin_dur])
    xticklabels(xticks*bin_dur-pre)
    title(bat_nms(i,:),'Color',bat_clr(i,:))

end

%% Subtract PSTH
load(fullfile(cd,'PSTH_table_bout_takeoff.mat'))
load(fullfile(cd,'PSTH_table_bout_takeoff_shuffle.mat'))
%=== FIGURE: PSTH of flight bouts for a bat vs any other bats
figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('Shuffled PSTH of bout (bin size = %.1fsec, smoothing window = %d bins, data=%s, method=%s)',bin_dur, smoothing, psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')
 
for i=1:n_tags
    idx = strcmp(tbl.RespBat,bat_nms(i,:)); % table idx to extract target data
    idx_shuffle = strcmp(tbl_shuffle.RespBat,bat_nms(i,:)); % table idx to extract target data
    binsize = bin_dur * 1000; % (msec)
    ntrials = sum(idx); % number of trials
    ntrials_shuffle = sum(idx_shuffle); % number of trials
    triallen = (pre+post)*Fs+1; % length of a trial (samples)
    

    flight_times_concatenated = cell2mat(tbl.FlightTimes(idx));
    flight_times_shuffle_concatenated = cell2mat(tbl_shuffle.FlightTimes(idx_shuffle));
    
    nexttile()
    ax = gca;
    ph = psth_subtract(flight_times_concatenated, flight_times_shuffle_concatenated,...
        binsize, Fs, ntrials, ntrials_shuffle, n_ite, triallen, smoothing, ax);
    ph.FaceColor = bat_clr(i,:);
    ph.EdgeColor = 'none';
    xticklabels(xticks/1000-pre)
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_summed_corrected_%s.png',psth_type)))

%=== FIGURE: PSTH of flight bouts for every bat pairs
figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,sprintf('Corrected PSTH of bouts (bin size = %.1fsec, smoothing window = %d bins, data=%s, method=%s)',bin_dur, smoothing,psth_type,seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Firing Rate (Hz)')

for i=1:n_tags
    for j=1:n_tags
        nexttile()
        if i~=j
            idx = strcmp(tbl.TrigBat,bat_nms(i,:)) & strcmp(tbl.RespBat,bat_nms(j,:)); % table idx to extract target data
            idx_shuffle = strcmp(tbl_shuffle.TrigBat,bat_nms(i,:)) & strcmp(tbl_shuffle.RespBat,bat_nms(j,:)); % table idx to extract target data
            binsize = bin_dur * 1000; % (msec)
            ntrials = sum(idx); % number of trials
            ntrials_shuffle = sum(idx_shuffle); % number of trials
            triallen = (pre+post)*Fs+1; % length of a trial (samples)
            
        
            flight_times_concatenated = cell2mat(tbl.FlightTimes(idx));
            flight_times_concatenated_shuffle = cell2mat(tbl_shuffle.FlightTimes(idx_shuffle));
            
            ax = gca;
            ph = psth_subtract(flight_times_concatenated, flight_times_concatenated_shuffle,...
                binsize, Fs, ntrials, ntrials_shuffle, n_ite, triallen, smoothing, ax);
            ph.FaceColor = bat_clr(j,:);
            ph.EdgeColor = 'none';
            xticklabels(xticks/1000-pre)
        end
    end
end
saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('PSTH_allPairs_corrected_%s.png',psth_type)))


%% Raster plot
% Load data
data_type = 'bout_takeoff';
% data_type = 'bout_landing';
% data_type = 'bflying_takeoff';
% data_type = 'bkanding_takeoff';
load(fullfile(cd,sprintf('PSTH_table_%s.mat',data_type)))

ruster = cell(1,n_tags);

%=== FIGURE
for i=1:n_tags
    figure(); set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.5 0.45]);
%     tl = tiledlayout(1,n_tags,'TileSpacing','tight');
    title(sprintf('Ruster plot (bat=%s, data=%s, method=%s)',bat_nms(i,:),psth_type,seg_alg))
    xlabel('Time (s)'); ylabel('Trial')

    idx = strcmp(tbl.RespBat,bat_nms(i,:)); % table idx to extract target data
    flight_times = tbl.FlightTimes(idx);
    ntrials = sum(idx); % number of trials

    hold on
    % For all trials...
    cnt = 1;
    for iTrial = 1:ntrials
        
        flights = flight_times{iTrial}'; % Get all spikes of respective trial (samples)
        xspikes = repmat(flights,2,1);         % Replicate array
        yspikes = nan(size(xspikes));       % NaN array
        
        if ~isempty(yspikes)
            yspikes(1,:) = cnt-1;                % Y-offset for raster plot
            yspikes(2,:) = cnt;
%             yspikes(1,:) = iTrial-1;                % Y-offset for raster plot
%             yspikes(2,:) = iTrial;
            plot(xspikes(:,1), yspikes(:,1), 'Color', 'k')

            ts = zeros(1,(pre+post)*Fs+1); ts(flights) = 1;
            ruster{1,i} = [ruster{1,i}; [iTrial ts]];
            cnt = cnt+1;
        end        
    end
    hold off
    xticks(0:Fs*10:(pre+post)*Fs)
    xticklabels(xticks/Fs-pre)

    saveas(gcf,fullfile(cd,'Figs','PSTH',sprintf('ruster_%s_%s.png',data_type,bat_nms(i,:))))
    pause(0.3)
end
save(fullfile(cd,sprintf('ruster_array_%s.mat',data_type)),'ruster')




%% PSTH for initial flights of flight mode (all sessions)
% Load data
load(fullfile(cd,'bouts_all.mat'))
seg_alg = 'GWM';
data_type = 'bout_takeoff';
bmode_all_now = bmode_all;
% bmode_all_now = bflying_all;

% [btakeoff,~] = findchanges(bmode_all_now,0,1);
[~,btakeoff] = findchanges(bmode_all_now,0,1);

% PSTH properties
pre = 30; post = 30; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 5; % window bin size of smooth filter

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < size(btakeoff,1)
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,[1:i-1,i+1:end],i) = sum(ts,1); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end
% ts_one = ts_one - ts_one_shuffle;
% ts_leader = sum(ts_one,2); % sum for followers
% ts_follower = sum(ts_one,3); % sum for leaders

ts_leader = sum(ts_one,2); % sum for followers
ts_follower = sum(ts_one,3); % sum for leaders

% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

%=== FIGURE: PSTH for leaders
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_leader(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('Triggered by %s',bat_nms(i,:)))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_leader_%ds_all_%s_1.png',post,seg_alg))) % save
% close gcf


%=== FIGURE: PSTH for follower
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_follower(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_follower_%ds_all_%s_1.png',post,seg_alg))) % save
% close gcf


%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
%     tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
%         if i==j
%             nexttile(); xticks([]);yticks([])
%         else
            nexttile();
%             bar(N(:,i,tr_cnt),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            bar(movmean(N(:,i,j),sm_wd),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
            xlim([0 (pre+post)/bin_dur])
            xticklabels(xticks*bin_dur-pre)
            
%             tr_cnt = tr_cnt + 1;
%         end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end

saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_mode_%ds_all_%s_2.png',post,seg_alg))) % save
% close gcf



%% Export PSTH for UMAP
seg_alg = 'GWM';

% Load data
data_type = 'bout_takeoff'; load(fullfile(cd,'bouts_all.mat')); bmode_all_now = bmode_all; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% data_type = 'bout_landing'; load(fullfile(cd,'bouts_all.mat')); bmode_all_now = bmode_all; [~,btakeoff] = findchanges(bmode_all_now,0,1);
% data_type = 'bflying_takeoff'; load(fullfile(cd,'Extracted_All.mat'),'bflying_all'); bmode_all_now = bflying_all; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% data_type = 'bflying_landing'; load(fullfile(cd,'Extracted_All.mat'),'bflying_all'); bmode_all_now = bflying_all; [~,btakeoff] = findchanges(bmode_all_now,0,1);


% PSTH properties
pre = 20; post = 20; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 10; % window bin size of smooth filter

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < T_all
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,[1:i-1,i+1:end],i) = sum(ts,1); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end

% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

% Check if PSTHs are correct
ts_follower = sum(ts_one,3); % sum for leaders
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_follower(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
    xticklabels(xticks*bin_dur-pre)
end

%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

% Convert to table data
tbl = cell(n_tags^2,size(N,1)+2);
for i=1:n_tags
    for j=1:n_tags
        tbl{(i-1)*6+j,1} = bat_nms(i,:);
        tbl{(i-1)*6+j,2} = bat_nms(j,:);
        tbl((i-1)*6+j,3:end) = num2cell(movmean(N(:,i,j),sm_wd));
    end
end
tbl = cell2table(tbl,'VariableNames',[{'trigger','responder'},arrayfun(@num2str, 1:size(N,1), 'UniformOutput', 0)]);
tbl(strcmp(tbl.trigger,tbl.responder),:) = [];

% save as csv
writetable(tbl,fullfile(cd,sprintf('%s_umap.csv',data_type)))

%% Export PSTH for UMAP (devide into 3 parts)
seg_alg = 'GWM';

% Load data
data_type = 'bout_takeoff'; load(fullfile(cd,'bouts_all.mat'));

x = [1 4; 5 8; 9 13];
tbl_combined = table();
for d = 1:size(x,1)
bmode_all_now = bmode_all(exp_d >= x(d,1) & exp_d <= x(d,2),:); [btakeoff,~] = findchanges(bmode_all_now,0,1);
% data_type = 'bout_landing'; load(fullfile(cd,'bouts_all.mat')); bmode_all_now = bmode_all; [~,btakeoff] = findchanges(bmode_all_now,0,1);
% data_type = 'bflying_takeoff'; load(fullfile(cd,'Extracted_All.mat'),'bflying_all'); bmode_all_now = bflying_all; [btakeoff,~] = findchanges(bmode_all_now,0,1);
% data_type = 'bflying_landing'; load(fullfile(cd,'Extracted_All.mat'),'bflying_all'); bmode_all_now = bflying_all; [~,btakeoff] = findchanges(bmode_all_now,0,1);


% PSTH properties
pre = 20; post = 20; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 10; % window bin size of smooth filter

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < size(btakeoff,1)
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,[1:i-1,i+1:end],i) = sum(ts,1); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end

% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

% % Check if PSTHs are correct
% ts_follower = sum(ts_one,3); % sum for leaders
% N = zeros(length(bin_edges)-1,n_tags);
% for i=1:n_tags
%     N(:,i) = accumarray(ts_g',ts_follower(:,i));
% end
% 
% figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
% tl = tiledlayout(1,n_tags,'TileSpacing','tight');
% title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
% xlabel(tl,'Time (s)'); ylabel(tl,'Count')
% for i=1:n_tags
%     nexttile();
% %     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
%     bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
%     plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
%     xlim([0 (pre+post)/bin_dur])
%     title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
%     xticklabels(xticks*bin_dur-pre)
% end

%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

% Convert to table data
tbl = cell(n_tags^2,size(N,1)+3);
for i=1:n_tags
    for j=1:n_tags
        tbl{(i-1)*6+j,1} = {d};
        tbl{(i-1)*6+j,2} = bat_nms(i,:);
        tbl{(i-1)*6+j,3} = bat_nms(j,:);
        tbl((i-1)*6+j,4:end) = num2cell(movmean(N(:,i,j),sm_wd));
    end
end
tbl = cell2table(tbl,'VariableNames',[{'sessions', 'trigger','responder'},arrayfun(@num2str, 1:size(N,1), 'UniformOutput', 0)]);
tbl(strcmp(tbl.trigger,tbl.responder),:) = [];

tbl_combined = [tbl_combined; tbl];
end

% save as csv
writetable(tbl_combined,fullfile(cd,sprintf('%s_umap.csv',data_type)))

%% Visualization of UMAP embodding processed in python
% Import processed psth from local directory
tbl = tbl_combined;
umap = readtable(fullfile(cd,'bout_takeoff_umap_processed.csv'),'NumHeaderLines',1);
umap(:,1) = [];

umap.Properties.VariableNames = {'var1','var2'};
umap.trigger = tbl.trigger;
umap.responder = tbl.responder;
umap.sessions = tbl.sessions;

%=== FIGURE: UMAP projection of PSTH for followers
figure(); hold on
for i=1:n_tags
    x = umap.var1(strcmp(umap.responder,bat_nms(i,:)));
    y = umap.var2(strcmp(umap.responder,bat_nms(i,:)));
    scatter(x,y,[],bat_clr(i,:),'filled','DisplayName',bat_nms(i,:))
end
hold off
lgd = legend('Location','best');
title(lgd,'Followers')
xlabel('var1'); ylabel('var2');
title('UMAP projection of PSTH')

% Import processed psth from local directory
tbl = tbl_combined;
umap = readtable(fullfile(cd,'bout_takeoff_umap_processed.csv'),'NumHeaderLines',1);
umap(:,1) = [];

umap.Properties.VariableNames = {'var1','var2'};
umap.trigger = tbl.trigger;
umap.responder = tbl.responder;
umap.sessions = tbl.sessions;
umap = umap(1:30,:);

%=== FIGURE: UMAP projection of PSTH for followers
figure(); hold on
for i=1:n_tags
    x = umap.var1(strcmp(umap.responder,bat_nms(i,:)));
    y = umap.var2(strcmp(umap.responder,bat_nms(i,:)));
    scatter(x,y,[],bat_clr(i,:),'filled','DisplayName',bat_nms(i,:))
end
hold off
lgd = legend('Location','best');
title(lgd,'Followers')
xlabel('var1'); ylabel('var2');
title('UMAP projection of PSTH (early)')

% Import processed psth from local directory
tbl = tbl_combined;
umap = readtable(fullfile(cd,'bout_takeoff_umap_processed.csv'),'NumHeaderLines',1);
umap(:,1) = [];

umap.Properties.VariableNames = {'var1','var2'};
umap.trigger = tbl.trigger;
umap.responder = tbl.responder;
umap.sessions = tbl.sessions;
umap = umap(31:60,:);

%=== FIGURE: UMAP projection of PSTH for followers
figure(); hold on
for i=1:n_tags
    x = umap.var1(strcmp(umap.responder,bat_nms(i,:)));
    y = umap.var2(strcmp(umap.responder,bat_nms(i,:)));
    scatter(x,y,[],bat_clr(i,:),'filled','DisplayName',bat_nms(i,:))
end
hold off
lgd = legend('Location','best');
title(lgd,'Followers')
xlabel('var1'); ylabel('var2');
title('UMAP projection of PSTH (middle)')


% Import processed psth from local directory
tbl = tbl_combined;
umap = readtable(fullfile(cd,'bout_takeoff_umap_processed.csv'),'NumHeaderLines',1);
umap(:,1) = [];

umap.Properties.VariableNames = {'var1','var2'};
umap.trigger = tbl.trigger;
umap.responder = tbl.responder;
umap.sessions = tbl.sessions;
umap = umap(61:90,:);

%=== FIGURE: UMAP projection of PSTH for followers
figure(); hold on
for i=1:n_tags
    x = umap.var1(strcmp(umap.responder,bat_nms(i,:)));
    y = umap.var2(strcmp(umap.responder,bat_nms(i,:)));
    scatter(x,y,[],bat_clr(i,:),'filled','DisplayName',bat_nms(i,:))
end
hold off
lgd = legend('Location','best');
title(lgd,'Followers')
xlabel('var1'); ylabel('var2');
title('UMAP projection of PSTH (late)')

%%
%=== FIGURE: PSTH of each isolet
idx = umap.var1 > -2;
umap.group(idx) = 1; umap.group(~idx) = 2; 

for gr = 1:2
    idx = umap.group == gr;
    trigbat = umap.trigger(idx);
    respbat = umap.responder(idx);
    n_points = sum(idx);

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
    title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg))
    xlabel(tl,'Time (s)'); ylabel(tl,'Count')

    for i=1:n_points
        b1 = trigbat(i); b2 = respbat(i);
        b1_idx = find(strcmp(bat_nms,b1));
        b2_idx = find(strcmp(bat_nms,b2));
        
        nexttile((b1_idx-1)*6+b2_idx);
        bar(movmean(N(:,b1_idx,b2_idx),sm_wd),'FaceColor',bat_clr(b2_idx,:),...
            'EdgeColor',bat_clr(b2_idx,:),'BarWidth',1); hold on
        plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
        xlim([0 (pre+post)/bin_dur])
        xticklabels(xticks*bin_dur-pre)
        title(sprintf('%s - %s',b1{:},b2{:}))
    end
end


%% Flight number plot
tbl = [];
%=== FIGURE: flight number plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.15 1 0.8]);
tl = tiledlayout(2, ceil(length(sessions)/2), 'TileSpacing', 'tight');
xlabel(tl,'Time (sec)'); ylabel(tl,'# of flights'); 
 
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying','n_tags','bat_clr','bat_nms','Fs')    
    [btakeoff,~] = findchanges(bflying,0,1);

    nexttile(); hold on    
    for i=1:n_tags
        bname = bat_nms(i,:);
        y = normalize(cumsum(btakeoff(:,i)),'range');
        x = 1:length(y);
        xq = 1:length(y)/100:length(y);
        yy = spline(x,y,xq);
       
        plot(yy,'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
        
        tbl = [tbl; [{bname},num2str(ss),arrayfun(@num2str, yy, 'UniformOutput', 0)]];
        
    end
    xticklabels(xticks/Fs)
    title(sprintf('Day%d',ss)); legend('Location','northwest')
    hold off
    
    pause(0.3)
end

% save as csv
tbl = cell2table(tbl);
writetable(tbl,fullfile(cd,'cumFlight_umap.csv'))

%% Visualization
% Import processed psth from local directory
umap = readtable(fullfile(cd,'cumFlight_umap_processed.csv'),'NumHeaderLines',1);
umap(:,1) = [];
umap.Properties.VariableNames = {'var1','var2'};
umap.bat = tbl.tbl1;
umap.day = tbl.tbl2;

%=== FIGURE: UMAP projection of PSTH for followers
figure(); hold on
for i=1:n_tags
    x = umap.var1(strcmp(umap.bat,bat_nms(i,:)));
    y = umap.var2(strcmp(umap.bat,bat_nms(i,:)));
    scatter(x,y,[],bat_clr(i,:),'filled','DisplayName',bat_nms(i,:))
end
hold off
lgd = legend('Location','best');
title(lgd,'Followers')
xlabel('var1'); ylabel('var2');
title('UMAP projection of cumulative flights')


%%
%=== FIGURE: Cumlative flight plot of each isolet
idx = umap.var1 > -2;
umap.group(idx) = 1; umap.group(~idx) = 2; 

for gr = 1:2
    idx = umap.group == gr;
    trigbat = umap.trigger(idx);
    respbat = umap.responder(idx);
    n_points = sum(idx);

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
    title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg_now))
    xlabel(tl,'Time (s)'); ylabel(tl,'Count')

    for i=1:n_points
        b1 = trigbat(i); b2 = respbat(i);
        b1_idx = find(strcmp(bat_nms,b1));
        b2_idx = find(strcmp(bat_nms,b2));
        
        nexttile((b1_idx-1)*6+b2_idx);
        bar(movmean(N(:,b1_idx,b2_idx),sm_wd),'FaceColor',bat_clr(b2_idx,:),...
            'EdgeColor',bat_clr(b2_idx,:),'BarWidth',1); hold on
        plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
        xlim([0 (pre+post)/bin_dur])
        xticklabels(xticks*bin_dur-pre)
        title(sprintf('%s - %s',b1{:},b2{:}))
    end
end


%% PSTH for initial flights of flight mode (all sessions)
num_ite = 100;

% PSTH properties
pre = 30; post = 30; % length around trigger (sec)
bin_dur = 0.5; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 10; % window bin size of smooth filter

% Compute PSTH
ts_one_shuffle = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

[btakeoff,~] = findchanges(bmode_all_now,0,1);


for ite = 1:num_ite
    if mod(ite,10) == 0
        disp(ite)
    end
    btakeoff_perm = zeros(size(btakeoff));
    
    for i=1:n_tags
        btakeoff_perm(:,i) = btakeoff(randperm(length(btakeoff(:,i))),i);
    end
    
    
for i=1:n_tags
    trig = find(btakeoff_perm(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff_perm;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < T_all
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one_shuffle(:,[1:i-1,i+1:end],i) = ts_one_shuffle(:,[1:i-1,i+1:end],i) + squeeze(sum(ts,1)); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end
end
ts_one_shuffle = ts_one_shuffle / num_ite;
ts_leader_shuffle = sum(ts_one_shuffle,2); % sum for followers
ts_follower_shuffle = sum(ts_one_shuffle,3); % sum for leaders


% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

%=== FIGURE: PSTH for leaders
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_leader_shuffle(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('Triggered by %s',bat_nms(i,:)))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_leader_%ds_shuffle_1.png',post))) % save
% close gcf


%=== FIGURE: PSTH for follower
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_follower_shuffle(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_follower_%ds_shuffle_2.png',post))) % save
% close gcf


%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one_shuffle(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
%     tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
%         if i==j
%             nexttile(); xticks([]);yticks([])
%         else
            nexttile();
%             bar(N(:,i,tr_cnt),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            bar(movmean(N(:,i,j),sm_wd),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
            xlim([0 (pre+post)/bin_dur])
            xticklabels(xticks*bin_dur-pre)
            
%             tr_cnt = tr_cnt + 1;
%         end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end

saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_bout_%ds_shuffle3.png',post))) % save
% close gcf

% save(fullfile(cd,sprintf('PSTH_shuffle_-%ds_%ds.mat',pre,post)),...
%     'ts_one_shuffle','ts_follower_shuffle','ts_leader_shuffle','pre','post','bin_dur','ts_t','bin_edges','sm_wd')


%% PSTH for initial flights of flight mode (all sessions)
num_ite = 1000;

% PSTH properties
pre = 100; post = 150; % length around trigger (sec)
bin_dur = 1; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 5; % window bin size of smooth filter

% Compute PSTH
ts_one_shuffle2 = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

[btakeoff,~] = findchanges(bmode_all_now,0,1);


for ite = 1:num_ite
    if mod(ite,10) == 0
        disp(ite)
    end
    btakeoff_perm = zeros(size(btakeoff));
    
    for i=1:n_tags
        btakeoff_perm(:,i) = btakeoff(randperm(length(btakeoff(:,i))),i);
    end
    
    
for i=1:n_tags
    trig = find(btakeoff_perm(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff_perm;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < T_all
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one_shuffle2(:,[1:i-1,i+1:end],i) = ts_one_shuffle2(:,[1:i-1,i+1:end],i) + squeeze(sum(ts,1)); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end
end
ts_one_shuffle2 = ts_one_shuffle2 / num_ite;
ts_leader_shuffle2 = sum(ts_one_shuffle2,2); % sum for followers
ts_follower_shuffle2 = sum(ts_one_shuffle2,3); % sum for leaders


% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

%=== FIGURE: PSTH for leaders
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_leader_shuffle2(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('Triggered by %s',bat_nms(i,:)))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_leader_%ds_shuffle_1.png',post))) % save
% close gcf


%=== FIGURE: PSTH for follower
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_follower_shuffle2(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_follower_%ds_shuffle_2.png',post))) % save
% close gcf


%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one_shuffle2(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
%     tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
%         if i==j
%             nexttile(); xticks([]);yticks([])
%         else
            nexttile();
%             bar(N(:,i,tr_cnt),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            bar(movmean(N(:,i,j),sm_wd),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
            xlim([0 (pre+post)/bin_dur])
            xticklabels(xticks*bin_dur-pre)
            
%             tr_cnt = tr_cnt + 1;
%         end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end

saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_bout_%ds_shuffle3.png',post))) % save
% close gcf

%% Subtract shuffled    
% load(fullfile(cd,sprintf('PSTH_shuffle_-%ds_%ds.mat',pre,post))) 
load(fullfile(cd,'bouts_all.mat'))
seg_alg_now = 'GWM';
bmode_all_now = bmode_all{strcmp(seg_alg,seg_alg_now)};

[btakeoff,~] = findchanges(bmode_all_now,0,1);

% PSTH properties
pre = 30; post = 30; % length around trigger (sec)
bin_dur = 1; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 5; % window bin size of smooth filter

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < T_all
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,[1:i-1,i+1:end],i) = sum(ts,1); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end
% ts_one = ts_one - ts_one_shuffle;
% ts_leader = sum(ts_one,2); % sum for followers
% ts_follower = sum(ts_one,3); % sum for leaders

ts_one = ts_one;
ts_leader = sum(ts_one,2); % sum for followers
ts_follower = sum(ts_one,3); % sum for leaders



% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

%=== FIGURE: PSTH for leaders
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_leader(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg_now))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('Triggered by %s',bat_nms(i,:)))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_leader_%ds_all_%s_1.png',post,seg_alg_now))) % save
% close gcf


%=== FIGURE: PSTH for follower
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_follower(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg_now))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile();
%     bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    bar(movmean(N(:,i),sm_wd),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1); hold on
    plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
    xlim([0 (pre+post)/bin_dur])
    title(sprintf('%s',bat_nms(i,:)),'FontSize',12,'Color',bat_clr(i,:))
    xticklabels(xticks*bin_dur-pre)
end
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_follower_%ds_all_%s_1.png',post,seg_alg_now))) % save
% close gcf


%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,sprintf('PSTH of takeoffs (method=%s)',seg_alg_now))
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
%     tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
%         if i==j
%             nexttile(); xticks([]);yticks([])
%         else
            nexttile();
%             bar(N(:,i,tr_cnt),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            bar(movmean(N(:,i,j),sm_wd),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
            xlim([0 (pre+post)/bin_dur])
            xticklabels(xticks*bin_dur-pre)
            
%             tr_cnt = tr_cnt + 1;
%         end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end

saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_mode_%ds_all_%s_2.png',post,seg_alg_now))) % save
% close gcf



%% Comparison of PSTH across sessions
% PSTH properties
pre = 50; post = 150; % length around trigger (sec)
bin_dur = 1; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 5; % window bin size of smooth filter


% psth_all = zeros(length(bin_edges)-1,n_tags-1,n_tags,length(sessions)); % time x trigger (leader) x bat x session
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder,'bouts.mat'))
    load(fullfile(sessions(ss).folder,sessions(ss).name),'T','bat_nms','Fs')

[btakeoff,~] = findchanges(bmode_now,0,1);

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags,n_tags); % time x follower (triggered) bats x leader (trigger) bat
% ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1); % trig x time x bat(n-1)
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 && trig(j)+post*Fs < T
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
%             ts(j,:,:) = btakeoff(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,[1:i-1,i+1:end],i) = sum(ts,1); % time x triggered bats x trigger bats: triggered counts for each bat
%     ts_all(:,i) = sum(sum(ts,1),3); % time x trigger bats: triggered counts for all bats
end

% Convert to binned data
ts_g = zeros(size(ts_t)); % bin id of each timepoint
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end


%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
%     tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
%         if i==j
%             nexttile(); xticks([]);yticks([])
%         else
            nexttile();
%             bar(N(:,i,tr_cnt),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            bar(movmean(N(:,i,j),sm_wd),'FaceColor',bat_clr(j,:),'EdgeColor',bat_clr(j,:),'BarWidth',1); hold on
            plot(repmat(pre/bin_dur,[1 2]),ylim,'k--'); hold off
            xlim([0 pre+post])
            xticklabels(xticks*bin_dur-pre)
            
%             tr_cnt = tr_cnt + 1;
%         end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end


% 
% for i=1:n_tags
%     psth_all(:,:,i,ss) = movmean(N(:,[1:i-1,i+1:end],i),sm_wd);
% end
end

%%
% PSTH properties
pre = 50; post = 150; % length around trigger (sec)
bin_dur = 1; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)
sm_wd = 5; % window bin size of smooth filter


% for i=1:n_tags % compute psth for i-the bat
for i=6:6
    ts_all = zeros((pre+post)*Fs+1,n_tags-1,length(sessions)); % psth for every pairs
    
    for ss=1:length(sessions)
        load(fullfile(sessions(ss).folder,'bouts.mat'))
        load(fullfile(sessions(ss).folder,sessions(ss).name),'T','bat_nms','Fs')
        [btakeoff,~] = findchanges(bmode_now,0,1);
        
        cnt = 1; % count
        for j=1:n_tags % set trigger of psth for every other bats
            
            if i~=j % exclude the pair of same bat
                trig = find(btakeoff(:,j)); % use j-th bat as trigger of psth
                bat_trigger = bat_nms(j,:); % name of j-th bat
                
                ts = zeros(length(trig),(pre+post)*Fs+1); % trigger x time
                for tr=1:length(trig)
                    if trig(tr)-pre*Fs > 0 && trig(tr)+post*Fs < T
                        ts(tr,:) = btakeoff(trig(tr)-pre*Fs:trig(tr)+post*Fs,i); % event for i-th bat around trigger by j-th bat
                    end
                end
                ts_all(:,cnt,ss) = sum(ts,1); % i-th bat psth around j-th bat flight for session ss
            
                cnt = cnt+1;
            end
        end
    end
    
    % Convert to binned data
    ts_g = zeros(size(ts_t)); % bin id of each timepoint
    for bb=1:length(bin_edges)-1
        if bb < length(bin_edges)-1
            ts_g(ts_t>=bin_edges(bb) & ts_t<bin_edges(bb+1)) = bb;
        else
            ts_g(ts_t>=bin_edges(bb) & ts_t<=bin_edges(bb+1)) = bb;
        end
    end
    
    N = zeros(length(bin_edges)-1,n_tags-1,length(sessions)); % binned psth. time x trigger bat x sessions
    for j=1:n_tags-1
        for ss=1:length(sessions)
            N(:,j,ss) = accumarray(ts_g',ts_all(:,j,ss));
        end
    end
    
    %=== FIGURE: PSTH for i-th bat bouts
    X = reshape(N,[size(N,1) size(N,2)*size(N,3)]);
    figure();
    imagesc(movmean(X,5)')
end
 
%%













%% Granger causality test of cumulative bout number
gctest_h = ones(n_tags,n_tags,length(sessions));
gctest_p = ones(n_tags,n_tags,length(sessions));
gctest_stat = ones(n_tags,n_tags,length(sessions));
gctest_c = ones(n_tags,n_tags,length(sessions));
gctest_kpss = ones(n_tags,n_tags,length(sessions));

for ss=1:length(sessions)
    tic
    load(fullfile(sessions(ss).folder,'bouts.mat'))
%     bmode = bmode_all(exp_d == ss,:); % flight mode = 1
    [btakeoff,~] = findchanges(bmode_now,0,1); % initial flights of flight mode
    
    max_iniflight = 1;
    for i=1:n_tags
        iniflight = find(btakeoff(:,i),1,'first');
        if max_iniflight < iniflight
            max_iniflight = iniflight; % set time zero to the maximum value of initial flights of first bauts
        end
    end

    for i=1:n_tags
        for j=1:n_tags
            if i~=j
                
%                 y1 = cumsum(btakeoff(max_iniflight:end,i)); % cause
%                 y2 = cumsum(btakeoff(max_iniflight:end,j)); % reponse
                y1 = double(btakeoff(max_iniflight:end,i)); % cause
                y2 = double(btakeoff(max_iniflight:end,j)); % reponse
                
                gctest_kpss(i,j,ss) = kpsstest(y1) & kpsstest(y2);

                [h,p,stat,cvalue] = gctest(y1,y2,NumLags=5000);
                gctest_h(i,j,ss) = h; gctest_p(i,j,ss) = p; gctest_stat(i,j,ss) = stat; gctest_c(i,j,ss) = cvalue;
            end
        end
    end
    toc
end

save(fullfile(cd,'gctest.mat'),'gctest_h','gctest_p','gctest_stat','gctest_c','gctest_kpss')

%%
y1 = cumsum(btakeoff(max_iniflight:end,i)); % cause
y2 = cumsum(btakeoff(max_iniflight:end,j)); % reponse

D1 = LagOp({1,-1},'Lags',[0,1]);
dy1 = filter(D1,y1);
dy2 = filter(D1,y2);

figure();
plot(y1); hold on; plot(y2); hold off

figure();
plot(dy1); hold on; plot(dy2); hold off

%% Cumulative flight mode plot for all sessions
cumflight_all = []; % cumulative flights for all sessions
cumflight_now = zeros(1,n_tags);
for ss=1:length(sessions)
    bmode_now = bmode_all_now(exp_d == ss,:); % flight mode = 1
    [btakeoff,~] = findchanges(bmode_now,0,1); % initial flights of flight mode
    
    max_iniflight = 1;
    for i=1:n_tags
        iniflight = find(btakeoff(:,i),1,'first');
        if max_iniflight < iniflight
            max_iniflight = iniflight; % set time zero to the maximum value of initial flights of first bauts
        end
    end

    % concatenate for all sessions
    cumflight_session = zeros(size(btakeoff(max_iniflight:end,:)));
    for i=1:n_tags
        cumflight_session(:,i) = cumsum(btakeoff(max_iniflight:end,i)); % cumulative number of flight mode
    end
    cumflight_all = [cumflight_all; cumflight_now+cumflight_session];
    cumflight_now = cumflight_all(end,:);
end

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.15 1 0.8]);
hold on
for i=1:n_tags
    plot(cumflight_all(:,i)/max(cumflight_all(:,i)),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
%     plot(cumflight_all(:,i),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
end
hold off
xticklabels(xticks/Fs)
title('Cumulative flight mode');
xlabel('Time (sec)')
ylabel('Normalized count')
legend('Location','northwest')

%% Distance of cumulative flight mode plot
dist_cumflight = zeros(size(cumflight_all,1),n_tags,n_tags);
for i=1:n_tags
    for j=1:n_tags
%         bat1 = bat_pairs(i,1); bat2 = bat_pairs(i,2);
        bat1 = i; bat2 = j;
%         dist_cumflight(:,i,j) = cumflight_all(:,bat1)./max(cumflight_all(:,bat1)) - cumflight_all(:,bat2)./max(cumflight_all(:,bat2));
    dist_cumflight(:,i,j) = cumflight_all(:,bat1) - cumflight_all(:,bat2);
    end
end

%=== FIGURE
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl = tiledlayout(n_tags, 1, 'TileSpacing', 'tight');
title(tl,'Distance of cumulative flight mode');
xlabel(tl,'Time (sec)'); ylabel(tl,'Distance'); 

% for i=1:size(bat_pairs,1)
for i=1:n_tags
    nexttile()
    hold on
    for j=1:n_tags
        if i~=j
            plot(dist_cumflight(:,i,j),'Color',bat_clr(j,:),'DisplayName',sprintf('%s - %s',bat_nms(i,:),bat_nms(j,:)))
%     plot(cumflight_all(:,i),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
        end
    end
    legend('Location','northwest')
    xticklabels(xticks/Fs)
    hold off
end

hold off




%% IFI for all flights of any bat
[btakeoff,~] = findchanges(bflying_all,0,1);
btakeoff_any = any(btakeoff,2);

bin_len = figopt.ifibin; % bin duration (timepoints)
bin_max = figopt.ifimax; % maximum ifi to count for IFI return plot
IFIedges = 0:bin_len:bin_max; 
IFI_any = diff(find(btakeoff_any));
IFI_any2 = IFI_any(IFI_any<bin_max);
ifi_now = IFI_any2(1:end-1);
ifi_next = IFI_any2(2:end);
[IFI_return_any,~,~] = histcounts2(ifi_now,ifi_next,IFIedges,IFIedges);

%=== FIGURE: IFI plot & IFI return plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.3 0.3 0.35 0.35]);
tl = tiledlayout(1, 2, 'TileSpacing', 'tight');
title(tl,'IFI plots for all flights (all sessions)')
% title(tl,'Inter-flight interval'); xlabel(tl,'Time (min)'); ylabel(tl,'Count')
hedge = figopt.hedges;
% IFI plot
nexttile(); histogram(IFI_any,hedge);
% title('IFI plot: all flights','FontSize',12)
xlabel('Time (sec)'); xticks(ifi_xticks); xticklabels(ifi_xticklbls); xlim([0,ifi_xlim])
ylabel('Count')
hold on; plot(repmat(median(IFI_any),[1,2]),ylim,'k--'); hold off
legend('',sprintf('Median: %.0fsec',median(IFI_any)/100));

%=== FIGURE: IFI return plot
nexttile(); imagesc(IFI_return_any);
set(gca,'YDir','normal')
% title('IFI return plot: all flights','FontSize',12)
xlabel('Time (sec)'); ylabel('Time (sec)')
xticks(figopt.ifiedges); xticklabels(figopt.ifilbls); 
yticks(figopt.ifiedges); yticklabels(figopt.ifilbls);
colormap jet; colorbar


%% Xcorr of takeoff timings
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name))
    [btakeoff,~] = findchanges(bflying,0,1);
    xcorr_btakeoff(btakeoff,n_tags,t,bat_nms,bat_clr,figopt)
    
    saveas(gcf,fullfile(cd,'Figs',['xcorr_takeoff_' num2str(ss) '.png']))
    close gcf
end

%% Xcorr of flight bouts
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'n_tags','t','bat_nms','bat_clr')
    load(fullfile(sessions(ss).folder, 'bouts.mat'))
    [btakeoff,~] = findchanges(bmode_now,0,1);
    xcorr_btakeoff(btakeoff,n_tags,t,bat_nms,bat_clr,figopt)
    
%     saveas(gcf,fullfile(cd,'Figs',['xcorr_takeoff_' num2str(ss) '.png']))
%     close gcf
end

%% Xcorr of takeoff timings (all sessions)
[btakeoff,~] = findchanges(bflying_all,0,1);
xcorr_btakeoff(btakeoff,n_tags,t,bat_nms,bat_clr,figopt)

%% Xcorr of flight bouts (all sessions)
% Parameters for xcorr
figopt.cc_bin = 0.5; % bin duration for cross-correlogram (sec)
figopt.cc_maxlag = 30; % maximum lag for cross-correlogram (sec)
figopt.cc_width = 4; % window size of mask for convlution (sec)
% figopt.cc_mask = ones(1,figopt.cc_width/figopt.cc_bin)/(figopt.cc_width/figopt.cc_bin); % mask for convolution smoothing
figopt.cc_mask = ones(1,figopt.cc_width/figopt.cc_bin); % mask for convolution smoothing
figopt.cc_ite = 100; % number of iteration for shift predictor

% Compute xcorr
load(fullfile(cd, 'Extracted_All.mat'),'n_tags','t_all','bat_nms','bat_clr')
load(fullfile(cd,'bouts_all.mat'),'bmode_all')
seg_alg_now = 'GWM';
bmode_all_now = bmode_all{strcmp(seg_alg,seg_alg_now)};
[btakeoff,~] = findchanges(bmode_all_now,0,1);
% xcorr_btakeoff(btakeoff,n_tags,t_all,bat_nms,bat_clr,figopt)

% options
maxlag = figopt.cc_maxlag; % maximum time to compute cc (sec)
mask = figopt.cc_mask;
bin_dur = figopt.cc_bin; % bin duration (sec)
bin_edges = 0:bin_dur:t_all(end);

% Compute xcorrs
xcorr_raw = zeros(2*maxlag/bin_dur+1,n_tags^2);

for i=1:n_tags^2
    id_col = i-n_tags*floor((i-1)/n_tags); % index of column of correlogram
    id_row = ceil(i/n_tags); % index of row of correlogram
    if id_col > id_row
        xcorr_raw(:,i) = zeros(size(xcorr_raw,1),1);
    
    else
        X = btakeoff(:,id_row); Y = btakeoff(:,id_col); % takeoffs of pairwise bats to compute xcorr
        N = zeros(length(bin_edges)-1,2);
        X_takeoff = t_all(X); Y_takeoff = t_all(Y);
        [N(:,1),~] = histcounts(X_takeoff,bin_edges);
        [N(:,2),~] = histcounts(Y_takeoff,bin_edges);
        [xcorr_raw(:,i), lags] = xcorr(N(:,1),N(:,2),maxlag/bin_dur);
    end
end

%=== FIGURE: Trajectories individual bats
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(n_tags, n_tags, 'TileSpacing', 'tight');
ylabel(tl,'Count'); xlabel(tl,'Time (sec)')
for i=1:n_tags^2
    nexttile(); 
    plot(bin_dur * lags,conv(xcorr_raw(:,i),mask,'same'),'k-'); hold on
    plot([0 0],ylim,'k--'); hold off
    xlim([-maxlag maxlag])
    if mod(i,n_tags)==1
        ylabel(bat_nms((i-1)/n_tags+1,:),'FontSize',12,'Color',bat_clr((i-1)/n_tags+1,:))
    end
    if i<=n_tags
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    end
    if i<=n_tags*(n_tags-1)
       xticks([]);
    end
end

%% Computing shift predictor
[btakeoff,~] = findchanges(bflying_all,0,1);
IFI_any = diff(find(any(btakeoff,2))); % IFIs of all flights

% options
bin_dur = figopt.cc_bin; % bin duration of xcorr (sec)
bin_edges = 0:bin_dur:t_all(end); % edges of xcorr (sec)
maxlag = figopt.cc_maxlag; % maximum time to compute cc (sec)
mask = figopt.cc_mask; % window size of smoothing
num_ite = figopt.cc_ite;


shft_dist = [2*prctile(IFI_any,95) T_all-2*prctile(IFI_any,95)]; % minimum and maximum shift


xcorr_shft = zeros(2*maxlag/bin_dur+1,n_tags^2);
for i=1:n_tags^2
    id_col = i-n_tags*floor((i-1)/n_tags); % index of column of correlogram
    id_row = ceil(i/n_tags); % index of row of correlogram
    if id_col > id_row
        xcorr_shft(:,i) = zeros(size(xcorr_shft,1),1);
    
    else
        X = btakeoff(:,id_row); Y = btakeoff(:,id_col); % takeoffs of pairwise bats to compute xcorr
        shft = randi(int32(shft_dist),[1 num_ite]);
        
        
        f = waitbar(0,'1','Name','Now computing...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(f,'canceling',0);
        msg = sprintf('Bat %d and bat %d',id_row,id_col);
        % iteration
        for j=1:num_ite
            % Check for clicked Cancel button
            if getappdata(f,'canceling'); break; end
            % Update waitbar and message
            waitbar(j/num_ite,f,msg);
            
            N = zeros(length(bin_edges)-1,2);
            Y = circshift(Y,shft);
            X_takeoff = t_all(X); Y_takeoff = t_all(Y);
            [N(:,1),~] = histcounts(X_takeoff,bin_edges);
            [N(:,2),~] = histcounts(Y_takeoff,bin_edges);
            [shft_pred, lags] = xcorr(N(:,1),N(:,2),maxlag/bin_dur);
            xcorr_shft(:,i) = xcorr_shft(:,i) + shft_pred;
        end
        delete(f)
    end
end

xcorr_shft = xcorr_shft/num_ite;

%=== FIGURE: plot shift predictor
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(n_tags, n_tags, 'TileSpacing', 'tight');
ylabel(tl,'Count'); xlabel(tl,'Time (sec)')

for i=1:size(xcorr_shft,2)
    nexttile(); plot(bin_dur * lags,conv(xcorr_shft(:,i),mask,'same'),'k-')
    if mod(i,n_tags)==1
        ylabel(bat_nms((i-1)/n_tags+1,:),'FontSize',12,'Color',bat_clr((i-1)/n_tags+1,:))
    end
    if i<=n_tags
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    end
    if i<=n_tags*(n_tags-1)
       xticks([]);
    end
end

%% Subtraction
xcorr_corrected = xcorr_raw - xcorr_shft;

%=== FIGURE: plot shift predictor
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(n_tags, n_tags, 'TileSpacing', 'tight');
ylabel(tl,'Count'); xlabel(tl,'Time (sec)')

for i=1:size(xcorr_corrected,2)
    nexttile(); plot(bin_dur * lags,conv(xcorr_corrected(:,i),mask,'same'),'k-')
    hold on; plot(xlim,zeros(size(xlim)),'k-'); hold off
    if mod(i,n_tags)==1
        ylabel(bat_nms((i-1)/n_tags+1,:),'FontSize',12,'Color',bat_clr((i-1)/n_tags+1,:))
    end
    if i<=n_tags
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    end
    if i<=n_tags*(n_tags-1)
       xticks([]);
    end
end

%% Statistical test
palpha = 0.05;
p_xcorr = zeros(size(xcorr_raw));
for i=1:size(xcorr_raw,2)
    id_col = i-n_tags*floor((i-1)/n_tags); % index of column of correlogram
    id_row = ceil(i/n_tags); % index of row of correlogram
    if id_col > id_row
        p_xcorr(:,i) = zeros(size(p_xcorr,1),1);
    else
        lambda = xcorr_shft(:,i);
        X = xcorr_shft(:,i);
        p_xcorr(:,i) = min(poisscdf(X,lambda),1-poisscdf(X,lambda));
    end
end

%=== FIGURE: Trajectories individual bats
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(n_tags, n_tags, 'TileSpacing', 'tight');
ylabel(tl,'Count'); xlabel(tl,'Time (sec)')
for i=1:n_tags^2
    nexttile(); plot(bin_dur * lags,conv(xcorr_raw(:,i),mask,'same'),'k-')
    if mod(i,n_tags)==1
        ylabel(bat_nms((i-1)/n_tags+1,:),'FontSize',12,'Color',bat_clr((i-1)/n_tags+1,:))
    end
    if i<=n_tags
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    end
    if i<=n_tags*(n_tags-1)
       xticks([]);
    end
end

%% Flight segmentation
for ss=1:length(sessions)
    flight_segmentation(fullfile(sessions(i).folder, sessions(i).name))
    saveas(gcf,fullfile(cd,'Figs',['flight_seg_day' num2str(ss) '.png']))
    
    pause(0.5)
    close gcf
end

%% Flight segmentation with w-xcorr
% for ss=1:length(sessions)
for ss=6:6
load(fullfile(sessions(ss).folder, sessions(ss).name),'n_tags','f_num','bat_nms','bflying','t')
load(fullfile(sessions(ss).folder, 'bpos.mat'),'bpos')

% Parameters for wxcorr
lagmax = 5; % (sec)
wsize = 5; % (sec)
winc = 1; % window increment (sec)

%=== FIGURE: Velocity and flight segmentation
n_ax = 1;

for i = 1:n_tags
    bnow = i;

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,1,'TileSpacing','tight');
    title(tl,sprintf('%s, %d flights',bat_nms(i,:),f_num(i)))

    ax(n_ax) = nexttile; n_ax = n_ax+1;
    % Flight plot
    area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    plot(t,bpos(:,i),'k--');
    ylabel('Position (cluster ID)'); hold off;
    legend({'Fly','Position'});

    % wxcorr plot
    for j=1:n_tags
        bpair = j;
        if strcmp(bat_nms(bnow,:),bat_nms(bpair,:))
        else
            [lag_time,twin,xcl] = timewindow_xcorr(bflying(:,bnow),bflying(:,bpair),Fs,wsize,winc,lagmax);
            
            ax(n_ax) = nexttile;
            imagesc(xcl'); 
            ylabel(sprintf('%s - %s',bat_nms(bnow,:),bat_nms(bpair,:)));
            xticks([]); ax(n_ax).TickLength = [0 0];
            yticks([1,lagmax*Fs+1,2*lagmax*Fs+1]); yticklabels((yticks-lagmax*Fs-1)/Fs);

            n_ax = n_ax+1;
        end
    end
    linkaxes(ax,'x');
    xlabel('Time(s)');
%     xticks(0:30*60/winc:T/(Fs*winc)); xticklabels(xticks/(60/winc));
    pause(0.5)
%     saveas(gcf,fullfile(cd,'Figs',sprintf('flight_wxcorr_day%d_b%d.png',ss,i)))
%     close gcf
end
end

%% Flight segmentation with w-xcorr (takeoff)
for ss=1:length(sessions)

load(fullfile(sessions(ss).folder, sessions(ss).name))
load(fullfile(sessions(ss).folder, 'bpos.mat'),'bpos')
[btakeoff,~] = findchanges(bflying,0,1);

% Parameters for wxcorr
lagmax = 30; % (sec)
wsize = 30; % (sec)
winc = 1; % window increment (sec)

%=== FIGURE: Velocity and flight segmentation
n_ax = 1;

for i = 1:n_tags
    bnow = i;

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,1,'TileSpacing','tight');
    title(tl,sprintf('%s, %d flights',bat_nms(i,:),f_num(i)))

    ax(n_ax) = nexttile; n_ax = n_ax+1;
    % Flight plot
    area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    plot(t,bpos(:,i),'k--');
    ylabel('Position (cluster ID)'); hold off;
    legend({'Fly','Position'});

    % wxcorr plot
    for j=1:n_tags
        bpair = j;
        if strcmp(bat_nms(bnow,:),bat_nms(bpair,:))
        else
            [lag_time,twin,xcl] = timewindow_xcorr(btakeoff(:,bnow),btakeoff(:,bpair),Fs,wsize,winc,lagmax);
            
            ax(n_ax) = nexttile;
            imagesc(xcl');            
            ylabel(sprintf('%s - %s',bat_nms(bnow,:),bat_nms(bpair,:)));
            xticks([]); ax(n_ax).TickLength = [0 0];
            yticks([1,lagmax*Fs+1,2*lagmax*Fs+1]); yticklabels((yticks-lagmax*Fs-1)/Fs);

            n_ax = n_ax+1;
        end
    end
    linkaxes(ax,'x');   xlabel('Time(s)');
%     xticks(0:30*60/winc:T/(Fs*winc)); xticklabels(xticks/(60/winc));
    pause(0.5)
    saveas(gcf,fullfile(cd,'Figs',sprintf('flight_wxcorr_tkoff_day%d_b%d.png',ss,i)))
    close gcf
end
end

%% Flight segmentation with w-xcorr (velocity)
for ss=2:length(sessions)

load(fullfile(sessions(ss).folder, sessions(ss).name))
load(fullfile(sessions(ss).folder, 'bpos.mat'),'bpos')

% Parameters for wxcorr
lagmax = 30; % (sec)
wsize = 30; % (sec)
winc = 1; % window increment (sec)

%=== FIGURE: Velocity and flight segmentation
n_ax = 1;

for i = 1:n_tags
    bnow = i;

    figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(n_tags,1,'TileSpacing','tight');
    title(tl,sprintf('%s, %d flights',bat_nms(i,:),f_num(i)))

    ax(n_ax) = nexttile; n_ax = n_ax+1;
    % Flight plot
    area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    plot(t,bpos(:,i),'k--');
    ylabel('Position (cluster ID)'); hold off;
    legend({'Fly','Position'});

    % wxcorr plot
    for j=1:n_tags
        bpair = j;
        if strcmp(bat_nms(bnow,:),bat_nms(bpair,:))
        else
            [lag_time,twin,xcl] = timewindow_xcorr(v_abs(:,bnow),v_abs(:,bpair),Fs,wsize,winc,lagmax);
            
            ax(n_ax) = nexttile;
            imagesc(xcl');            
            ylabel(sprintf('%s - %s',bat_nms(bnow,:),bat_nms(bpair,:)));
            xticks([]); ax(n_ax).TickLength = [0 0];
            yticks([1,lagmax*Fs+1,2*lagmax*Fs+1]); yticklabels((yticks-lagmax*Fs-1)/Fs);

            n_ax = n_ax+1;
        end
    end
    linkaxes(ax,'x');   xlabel('Time(s)');
%     xticks(0:30*60/winc:T/(Fs*winc)); xticklabels(xticks/(60/winc));
    pause(0.5)
    saveas(gcf,fullfile(cd,'Figs',sprintf('flight_wxcorr_vel_day%d_b%d.png',ss,i)))
    close gcf
end
end


%% Take-offs PSTH triggered by one bat
for ss = 1:length(sessions)
load(fullfile(sessions(ss).folder, sessions(ss).name))
[btakeoff,~] = findchanges(bflying,0,1);

% PSTH properties
pre = 100; post = 300; % length around trigger (sec)
bin_dur = 5; % bin duration (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags-1,n_tags); % time x triggered bats x trigger bat
ts_all = zeros((pre+post)*Fs+1,n_tags); % time x trigger bat

for i=1:n_tags
% for i=1:1
    trig = find(btakeoff(:,i));
    ts = zeros(length(trig),(pre+post)*Fs+1,n_tags-1);
    btakeoff_others = btakeoff;
    btakeoff_others(:,i) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 & trig(j)+post*Fs < T
%             sum(btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:))
            ts(j,:,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs,:);
        end
    end
    
    ts_one(:,:,i) = sum(ts,1); % triggered counts for each bat
    ts_all(:,i) = sum(sum(ts,3),1); % triggered counts for all bats
end    

% Convert to binned data
ts_g = zeros(size(ts_t));
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end

%=== FIGURE: PSTH for all bats' takeoffs
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_all(:,i));
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile(); bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1);
    title(sprintf('Triggered by %s',bat_nms(i,:)))
    xticklabels(xticks*bin_dur-pre)
end

%===SAVE
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_takeoff_%ds_day%d_1.png',pre+post,ss)))
close gcf

%=== FIGURE: PSTH for all pairs
N = zeros(length(bin_edges)-1,n_tags,n_tags-1);
for i=1:n_tags
    for j=1:n_tags-1
        N(:,i,j) = accumarray(ts_g',ts_one(:,j,i)); % 1st col: binned time, 2nd col: trigger, 3rd col: triggered bats
    end
end

figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(n_tags,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags % trigger
    tr_cnt = 1; % Count for triggered bats
    for j=1:n_tags % triggered bat
        if i==j
            nexttile(); xticks([]);yticks([])
        else
            nexttile(); bar(N(:,i,tr_cnt),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1);
            xticklabels(xticks*bin_dur-pre)
            
            tr_cnt = tr_cnt + 1;
        end

        if j==1
            ylabel(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        end
        if i==1
            title(bat_nms(j,:),'FontSize',12,'Color',bat_clr(j,:))
        end
    end
end

%===SAVE
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_takeoff_%ds_day%d_2.png',pre+post,ss)))
close gcf

end

%% Take-offs PSTH triggered by all other bats
for ss = 1:length(sessions)

load(fullfile(sessions(ss).folder, sessions(ss).name))
[btakeoff,~] = findchanges(bflying,0,1);

% PSTH properties
pre = 100; post = 300; % length around trigger (sec)
ts_t = -1*pre:1/Fs:post; % time around trigger (sec)
bin_dur = 5; % bin duration (sec)
bin_edges = ts_t(1):bin_dur:ts_t(end); % edges of bin (sec)

% Compute PSTH
ts_one = zeros((pre+post)*Fs+1,n_tags); % time x triggered bats

for i=1:n_tags
% for i=1:1
    trig = find(any(btakeoff(:,[1:i-1,i+1:end]),2)); % Triggered by takeoffs of any other bats
    ts = zeros(length(trig),(pre+post)*Fs+1); % trig x time x triggered bats
    btakeoff_others = btakeoff;
    btakeoff_others(:,[1:i-1,i+1:end]) = []; % remove the column for the trigger bat
    for j=1:length(trig)
        if trig(j)-pre*Fs > 0 & trig(j)+post*Fs < T
            ts(j,:) = btakeoff_others(trig(j)-pre*Fs:trig(j)+post*Fs);
        end
    end
    
    ts_one(:,i) = sum(ts,1); % triggered counts for each bat
end    

% Convert to binned data
ts_g = zeros(size(ts_t));
for i=1:length(bin_edges)-1
    if i < length(bin_edges)-1
        ts_g(ts_t>=bin_edges(i) & ts_t<bin_edges(i+1)) = i;
    else
        ts_g(ts_t>=bin_edges(i) & ts_t<=bin_edges(i+1)) = i;
    end
end
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    N(:,i) = accumarray(ts_g',ts_one(:,i));
end

%=== FIGURE: PSTH for all bats' takeoffs
figure(); set(gcf, 'units','normalized','outerposition',[0 0.3 1 0.4]);
tl = tiledlayout(1,n_tags,'TileSpacing','tight');
title(tl,'PSTH of takeoffs triggered by all other bats')
xlabel(tl,'Time (s)'); ylabel(tl,'Count')

for i=1:n_tags
    nexttile(); bar(N(:,i),'FaceColor',bat_clr(i,:),'EdgeColor',bat_clr(i,:),'BarWidth',1);
    title(bat_nms(i,:))
    xticklabels(xticks*bin_dur-pre)
end

%===SAVE
saveas(gcf,fullfile(cd,'Figs',sprintf('PSTH_allother_%ds_day%d.png',pre+post,ss)))
close gcf

end

%% Takeoffs triggered by reward
ss=1;
load(fullfile(sessions(ss).folder, sessions(ss).name))
load(fullfile(TTLs(ss).folder, TTLs(ss).name),'AnalogSignals')
rwds = AnalogSignals(:,1);


%% Window cross correlation
lagmax = 30; % (sec)
wsize = 30; % (sec)
winc = 1; % window increment (sec)
dd = 6; % day

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(size(bat_pairs,1), 1, 'TileSpacing', 'none');
title(tl,sprintf('Windowed cross-correlation (day%d): lag=%dsec, window=%dsec',dd,lagmax,wsize))
xlabel(tl,'Time (min)');

for i = 1:size(bat_pairs,1)
    b1 = bat_pairs(i,1); b2 = bat_pairs(i,2);
    % [lag_time,twin,xcl] = timewindow_xcorr(bflying_low(:,1),bflying_low(:,2),dwn_fs,20,1,5);
    [lag_time,twin,xcl] = timewindow_xcorr(bflying(:,b1),bflying(:,b2),Fs,wsize,winc,lagmax);

    ax(i) = nexttile; imagesc(xcl')
    ylabel(bat_pair_nms(i,:));
%     yticks([1,lagmax*Fs+1,2*lagmax*Fs+1]); yticklabels((yticks-lagmax*Fs-1)/Fs);
    yticks([])
    if i~=size(bat_pairs,1)
        xticks([])
    else
        xticks(0:30*60/winc:T/(Fs*winc)); xticklabels(xticks/(60/winc));
    end
    ax(i).TickLength = [0 0];
end
linkaxes(ax,'x');

%% Window xcorr and flying timings


%% Window cross correlation for position
maxlag = 1000;
wsize = 10;
wmask = [zeros(1,maxlag-wsize),ones(1,2*wsize+1),zeros(1,maxlag-wsize)];
winc = 100;
b1 = 2; b2 = 1;

if wsize < maxlag
    WC = zeros(floor(T/winc-2*maxlag),2*maxlag+1);
    for i=1:size(WC,1)
        idx = 1+(i-1)*winc;
        X =  bflying(idx:idx+2*maxlag,b1); X = X.*wmask';
        Y =  bflying(idx:idx+2*maxlag,b2);
        [wcc, lags] = xcorr(X,Y,maxlag); 
        WC(i,:) = wcc;
    end
else
    WC = zeros(floor(T/winc-wsize+1),2*maxlag+1);
    for i=1:size(WC,1)
        idx = 1+(i-1)*winc;
        X =  wcc_data(idx:idx+wsize-1,1,b1);
        Y =  wcc_data(idx:idx+wsize-1,1,b2);
        [wcc, lags] = xcorr(X,Y,maxlag); 
        WC(i,:) = wcc;
    end
end

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.5 1 0.35]);
imagesc(WC')
yticks(linspace(1,2*maxlag+1,2*maxlag/500+1)); yticklabels((yticks-maxlag-1)/500)

%%
figopt.wcc_wsize = 0.1; % window size (sec)
figopt.wcc_winc = 1; % window increment (sec)
figopt.wcc_maxlag = 1000; % maximum lag (csec)
figopt.wcc_laginc = 50; % lag increment (csec)

figopt.wcc_bin = 0.50; % bin duration for cross-correlogram (sec)
figopt.wcc_maxlag = 5; % maximum lag for cross-correlogram (sec)
figopt.wcc_wsize = 10;
figopt.ccwidth = 0.5; % window size of mask for convlution (sec)

lagmax = figopt.wcc_maxlag;
bin_dur = figopt.wcc_bin; % bin duration (sec)
wsize = figopt.wcc_wsize; % window size (sec), must be greater than lagmax

[btakeoff,~] = findchanges(bflying_all,0,1);
WC = zeros(floor((t_all(end)-wsize+1)/bin_dur),2*lagmax/bin_dur+1);

% convert to binned data (downsampling)
bin_edges = 0:bin_dur:t_all(end);
N = zeros(length(bin_edges)-1,n_tags);
for i=1:n_tags
    t_takeoff = t_all(btakeoff(:,i));
    [N(:,i),~] = histcounts(t_takeoff,bin_edges);
end

% Calculate cross-correlations for time i
for i=1:size(WC,1)
    X = N(i:i+(wsize-1)/bin_dur,1);
    Y = N(i:i+(wsize-1)/bin_dur,1);
    [wcc_takeoff, lags] = xcorr(X,Y,maxlag/bin_dur); 
    WC(i,:) = wcc_takeoff;
end


maxlag = figopt.ccmaxlag; % maximum time to compute cc (sec)
[xcorr_raw, lags] = xcorr(N,maxlag/bin_dur); 
mask = figopt.ccmask;


figopt.wcc_wsize = 10; % window size (csec)
figopt.wcc_winc = 10; % window increment (sec)
figopt.wcc_maxlag = 1000; % maximum lag (csec)
figopt.wcc_laginc = 50; % lag increment (csec)

wsize = figopt.wcc_wsize;
winc = figopt.wcc_winc;
lagmax = figopt.wcc_maxlag;
laginc = figopt.wcc_laginc;

WC = zeros(floor((T_all-2*lagmax-wsize+winc)/winc),2*lagmax/laginc+1); % windowed cross correlation
textprogressbar('Computing window cross correlation...   ')
for i=1:size(WC,1)
    textprogressbar(100*i/size(WC,1))
    idx = lagmax+1+(i-1)*winc;
    X = [zeros(lagmax,1); btakeoff(idx:idx+wsize-1,1); zeros(lagmax,1)];
    Y = btakeoff(idx-lagmax:idx+lagmax+wsize-1,2);
    wcc = xcorr(X,Y,lagmax);
    wcc = wcc(linspace(1,2*lagmax+1,2*lagmax/laginc+1));
    WC(i,:) = wcc;
end
textprogressbar('done')






%% Inter flight intervals
mkdir(fullfile(cd,'Figs','IFI'))

scale = {'log','linear'};

% Compute IFI
IFI_med = zeros(length(sessions),n_tags);
IFI_uppq = zeros(length(sessions),n_tags);
IFI_lowq = zeros(length(sessions),n_tags);
IFI_ret_xmean = zeros(length(sessions),n_tags);
IFI_ret_ymean = zeros(length(sessions),n_tags);
IFI_ret_xstd = zeros(length(sessions),n_tags);
IFI_ret_ystd = zeros(length(sessions),n_tags);

for i=1:2
    for ss=1:length(sessions)
    % for ss = 1:1
        load(fullfile(sessions(ss).folder, sessions(ss).name))
        [btakeoff,~] = findchanges(bflying,0,1);

        % IFI plot
        [IFI_med(ss,:), IFI_uppq(ss,:),IFI_lowq(ss,:)] = compute_IFI(btakeoff,bat_nms,bat_clr,Fs,figopt,ss,scale{i});
        saveas(gcf,fullfile(cd,'Figs','IFI',sprintf('IFI_%s_day%d.png',scale{i},ss)))

        % IFI_return plot
        [IFI_ret_xmean(ss,:), IFI_ret_ymean(ss,:), IFI_ret_xstd(ss,:), IFI_ret_ystd(ss,:)]...
            = compute_IFI_return(btakeoff,bat_nms, bat_clr,Fs,figopt,ss,scale{i});
        saveas(gcf,fullfile(cd,'Figs','IFI',sprintf('IFIreturn_%s_day%d.png',scale{i},ss)))

        close all
        pause(0.3)
    end
end


% IFI median plot across sessions
figure(); hold on
for i=1:n_tags
    ax = plot(IFI_med(:,i),'o-');
    set(ax,'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
end
hold off
legend(); xlim([1 length(sessions)])
xlabel('session'); ylabel('IFI median (sec)');
saveas(gcf,fullfile(cd,'Figs','IFI','IFI_med.png'))

% IFI return mean scatter
figure(); hold on
for i = 1:n_tags
    ax = scatter(IFI_ret_xmean(:,i),IFI_ret_ymean(:,i),'filled','o');
    set(ax,'MarkerEdgeColor',bat_clr(i,:),'DisplayName',bat_nms(i,:))
end
hold off
legend()
xlabel('IFI(n) mean (sec)'); ylabel('IFI(n+1) mean (sec)');
saveas(gcf,fullfile(cd,'Figs','IFI','IFI_ret_mean.png'))

% IFI return std scatter
figure(); hold on
for i = 1:n_tags
    ax = scatter(IFI_ret_xstd(:,i),IFI_ret_ystd(:,i),'filled','o');
    set(ax,'MarkerEdgeColor',bat_clr(i,:),'DisplayName',bat_nms(i,:))
end
hold off
legend()
xlabel('IFI(n) std (sec)'); ylabel('IFI(n+1) std (sec)');
saveas(gcf,fullfile(cd,'Figs','IFI','IFI_ret_std.png'))

close all

%% IFI for all sessions
scale = {'log','linear'};

[btakeoff,~] = findchanges(bflying_all,0,1);
for i=1:2

% IFI plot
[~,~,~] = compute_IFI(btakeoff,bat_nms,bat_clr,Fs,figopt,0,scale{i});
saveas(gcf,fullfile(cd,'Figs','IFI',sprintf('IFI_all_%s.png',scale{i})))

% IFI_return plot
[~, ~, ~, ~] = compute_IFI_return(btakeoff,bat_nms, bat_clr,Fs,figopt,0,scale{i});
saveas(gcf,fullfile(cd,'Figs','IFI',sprintf('IFIreturn_all_%s.png',scale{i})))

end

%% Flight number plot
%=== FIGURE: flight number plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.15 1 0.8]);
tl = tiledlayout(2, ceil(length(sessions)/2), 'TileSpacing', 'tight');
 xlabel(tl,'Time (sec)'); ylabel(tl,'# of flights'); 
 
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder, sessions(ss).name),'bflying','n_tags','bat_clr','bat_nms','Fs')    
    [btakeoff,~] = findchanges(bflying,0,1);

    nexttile(); hold on    
    for i=1:n_tags
%         plot(cumsum(btakeoff(:,i)),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
        plot(normalize(cumsum(btakeoff(:,i)),'range'),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
    end
    xticklabels(xticks/Fs)
    title(sprintf('Day%d',ss)); legend('Location','northwest')
    hold off
    
    pause(0.3)
end

%% Fligth bouts number plot
%=== FIGURE: flight bout number plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.15 1 0.8]);
tl = tiledlayout(2, ceil(length(sessions)/2), 'TileSpacing', 'tight');
xlabel(tl,'Time (sec)'); ylabel(tl,'# of flights'); 
for ss=1:length(sessions)
    load(fullfile(sessions(ss).folder,'bouts.mat')) % load bat flight bouts
    load(fullfile(sessions(ss).folder,sessions(ss).name),'bat_clr','bat_nms','Fs')

    [btakeoff,~] = findchanges(bmode_now,0,1); % initial flights of flight mode
    
    max_iniflight = 1;
    for i=1:n_tags
        iniflight = find(btakeoff(:,i),1,'first');
        if max_iniflight < iniflight
            max_iniflight = iniflight; % set time zero to the maximum value of initial flights of first bauts
        end
    end

    nexttile(); hold on
    for i=1:n_tags
        cumflight = cumsum(btakeoff(max_iniflight:end,i)); % cumulative number of flight mode
        plot(normalize(cumflight,'range'),'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
%         plot(cumflight,'Color',bat_clr(i,:),'DisplayName',bat_nms(i,:))
    end
    xticklabels(xticks/Fs)
    title(sprintf('Day%d',ss)); legend('Location','northwest')
    hold off
    
    pause(0.3)
end

%% Flights of other during pre/during/post bout



end