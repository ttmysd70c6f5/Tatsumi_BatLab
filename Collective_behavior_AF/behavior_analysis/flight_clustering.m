function flight_clustering()
% Hierarchical clustering of trajectories for all sessions

sessionFile = dir(fullfile(pwd,'*/*cdp*'));

%% Parameter
n_smp = 10; % number of downsampled points
n_dim = 3; % number of dimension (x,y,z)
cutoffs = [0.9:0.1:1.2]; % Cutoff value for clustering
min_trj = 5; % Minimum amounts of trajectory in a cluster

%% Interpolate trajectories
r_interp = [];
flights_all = [];
for ss = 1:length(sessionFile)
    % Load data
    sessionDir = sessionFile(ss).folder;
    BHV1_file = dir(fullfile(sessionDir, 'Ext_Beh*','Extracted_Behavior*'));                              % Preprocessed Behavioral Data
    load(fullfile(BHV1_file.folder,BHV1_file.name),'r');        % Processed Behavioral Data

    BHV3_file = dir(fullfile(sessionDir,'Analysis_TY','behaviors_extracted.mat'));
    load(fullfile(BHV3_file.folder,BHV3_file.name),'flights');

    % Concatenate
    r_interp_now = zeros(size(flights,1),n_dim*n_smp); % interpolated trajectory
    
    for nf = 1:size(flights,1)
        nb = flights.bat_id(nf); % Bat id
        r_flight = r(flights.ts_toff(nf):flights.ts_land(nf),:,nb); % raw trajectory
        sp_flight = Rate(flights.ts_toff(nf):flights.ts_land(nf),:); % firing rate
        fl_len_ts = flights.ts_land(nf) - flights.ts_toff(nf); % length of flight (timestamps)
        w_interp = 1 + fl_len_ts .* (0.5 + 0.5 .* (linspace(-n_smp/2,n_smp/2,n_smp)).^3 / (n_smp/2)^3); % weight for interpolation: sigmoid
        r_interp_now(nf,:) = reshape(interp1(r_flight, w_interp),[1,n_smp*3]);
    end

    r_interp = [r_interp; r_interp_now];
    flights_all = [flights_all;flights];
end

%% hierarchical clustering of trajecoties


[best_cutoff, num_good_clus] = optimize_prm_hc(r_interp, cutoffs, min_trj);

% cutoff = 1.1 resutled in the largest number of clusters
% cutoff = 1.2 resulted in only one cluster

%% Clustering with the best cutoff
% Parameter
n_smp = 10; % number of downsampled points
n_dim = 3; % number of dimension (x,y,z)
cutoff = 1.1; % Cutoff value for clustering
min_trj = 5; % Minimum amounts of trajectory in a cluster

% clustering
hc = struct([]);
% hc_clus = cell(n_tags,1);
% hc_good = cell(n_tags,1);
% hc_cts = cell(n_tags,1);
% hc_id = cell(n_tags,1);
% hc_good = cell(n_tags,1);

fig_cnt = 1;
for nb=1:n_tags
    % hierarchical clustering
    hc_tree = linkage(r_interp(flights_all.bat_id==nb,:),'single','euclidean'); % hierarchical cluster tree
    D = pdist(r_interp(flights_all.bat_id==nb,:));
    leafOrder = optimalleaforder(hc_tree,D);
    hc(nb).hc_clus = cluster(hc_tree,'cutoff',cutoff,'Criterion','inconsistent');
    [hc(nb).hc_cts,hc(nb).hc_id] = groupcounts(hc(nb).hc_clus);
    hc(nb).hc_good = hc(nb).hc_id(hc(nb).hc_cts > min_trj);

    fprintf('\n Number of clusters for bat %d: %d \n',nb,length(unique(hc(nb).hc_id)))
    fprintf('Number of clusters with more than %d flights for bat %d: %d \n',min_trj,nb,length(unique(hc(nb).hc_good)))

    % takeoff and landing position for good clusters
    for nclu = 1:length(hc(nb).hc_good)
        r_toff_ext = flights_all.r_toff(flights_all.bat_id==nb,:);
        r_land_ext = flights_all.r_land(flights_all.bat_id==nb,:);

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
    tl = tiledlayout(6,11,'TileSpacing','tight');

    flights_ext = flights_all(flights_all.bat_id==nb,:);
    r_interp_ext = r_interp(flights_all.bat_id==nb,:);
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
        title(tl,sprintf('Hierarchical clustering of trajectory: bat %d, cutoff=%0.1f, minimum number of trajectory=%d',nb,round(cutoff,1),min_trj))
        xlabel(tl,'x'); ylabel(tl,'y')
    end
    exportgraphics(f,fullfile('G:\My Drive\UCBerkeley\Data Analysis\Progress_Meeting\Progress_2023-04-12_Michael',...
            sprintf('Trajectory_allsession_bat%.png',nb)))
    fig_cnt = fig_cnt + 1;
    pause(0.1)
    close gcf
end


%% Visualize clustered trajectories
% for nb = 1:n_tags
%     %===FIGURE: visualize clustered trajectories for each bat
%     figure
%     set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%     tl = tiledlayout(6,11,'TileSpacing','tight');
% 
%     flights_ext = flights_all(flights_all.bat_id==nb,:);
%     r_interp_ext = r_interp(flights_all.bat_id==nb,:);
%     r_interp_ext = reshape(r_interp_ext,[size(r_interp_ext,1),n_smp,n_dim]);
%     for nc = 1:length(hc(nb).hc_good)
%         nexttile; hold on
% 
%         flights_clu_now = flights_ext(hc(nb).hc_clus == hc(nb).hc_good(nc),:);
%         r_interp_clu_now = r_interp_ext(hc(nb).hc_clus == hc(nb).hc_good(nc),:,:);
%         
% %         x = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
% %         y = zeros(size(flights_clu_now,1),length(1:.1:n_smp));
%         c = 1:.1:n_smp; c = c/max(c); c = 100.^c;
%         hold on
%         for nf = 1:size(flights_clu_now,1)
%             x = spline(1:n_smp,r_interp_clu_now(nf,:,1),1:.1:n_smp);
%             y = spline(1:n_smp,r_interp_clu_now(nf,:,2),1:.1:n_smp);
%             scatter(x,y,10,c,'filled','MarkerFaceAlpha',0.5)       
%         end
%         hold off
%         colormap parula
%         xlim(r_lim(1,:)); ylim(r_lim(2,:));
%         title(tl,sprintf('Hierarchical clustering of trajectory: bat %d, cutoff=%d, minimum number of trajectory=%d',nb,cutoff,min_trj))
%         xlabel(tl,'x'); ylabel(tl,'y')
%     end
% end

%%

figure(); dendrogram(hc_tree,150,'Reorder',leafOrder,'ColorThreshold','default')

%% 
figure
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl = tiledlayout(7,12,'TileSpacing','tight');

for nclu = 1:length(hc_good)
    
    flights_clu_now = flights_all(hc_clus == hc_good(nclu),:);
    
    nexttile
    hold on
    for nf = 1:size(flights_clu_now,1)
        nb = flights_clu_now.bat_id(nf);
        ts_toff = flights_clu_now.ts_toff(nf);
        ts_land = flights_clu_now.ts_land(nf);
        
        sessionDir = fullfile(pwd,flights_clu_now.rec_date(nf,:));
        
        % Load data
        BHV1_file = dir(fullfile(sessionDir, 'Ext_Beh*','Extracted_Behavior*'));                              % Preprocessed Behavioral Data
        load(fullfile(BHV1_file.folder,BHV1_file.name),'r','r_lim');        % Processed Behavioral Data
        
        x = r(ts_toff:ts_land,1,nb);
        y = r(ts_toff:ts_land,2,nb);
        scatter(x(:),y(:),5,'k','filled','MarkerFaceAlpha',0.15)
        xlim(r_lim(1,:)); ylim(r_lim(2,:));
    end
    hold off
end

%% 
% for nc = 1:length(hc_good)
%     flights_clu_now = flights_all(hc_clus == hc_good(nc),:);