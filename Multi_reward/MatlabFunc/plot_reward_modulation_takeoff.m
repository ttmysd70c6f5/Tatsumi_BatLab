function plot_reward_modulation_takeoff(bname_ephys,trjName,Flights,TrjCluster,spikeData,mtData,Layers)
% Specify the ephys bat
bb = find(strcmp(bname_ephys,mtData.bname));

% Load behavioral and neural features
n_bats = length(mtData.bname); % number of bats
oth_ids = [1:bb-1,bb+1:n_bats]; % index for other conspecifics
t_flight_sec = Flights(bb).t_flight_sec;
is_fd_land = Flights(bb).is_fd_land;
% closest_feeder = Flights(bb).closest_feeder(:,2);
delivered_reward = Flights(bb).reward;
n_land = length(t_flight_sec); % number of landings
cluster_idx = TrjCluster(bb).clu_idx; % flight cluster index
cluster_list = TrjCluster(bb).clu_list; % unique clusters
cluster_count = TrjCluster(bb).clu_sz; % number of flights in each cluster
n_cluster = length(cluster_count); % number of clusters
n_unit = spikeData.numGoodUnits; % number of single units
unit_ids = spikeData.good_units.cluster_id;

t_margin = [-1 3]; % time window to extract spike times around landings (sec)

for clu = 1:n_cluster
    n_trj = cluster_count(clu); % number of trajectory for the current cluster
        
    if n_trj >= 3 && cluster_list(clu) ~= -1 % Plot only the cluster with more than 2 trajectories
        idx_r1 = cluster_idx == cluster_list(clu) & is_fd_land & delivered_reward == 1; % Extract flights for reward 1
        idx_r2 = cluster_idx == cluster_list(clu) & is_fd_land & delivered_reward == 2; % Extract flights for reward 2
        if sum(any([idx_r1,idx_r2],2)) ~= 0
            t_onset_sec_r1 = t_flight_sec(idx_r1,1);
            t_end_sec_r1 = t_flight_sec(idx_r1,1);
        
            t_onset_sec_r2 = t_flight_sec(idx_r2,1);
            t_end_sec_r2 = t_flight_sec(idx_r2,1);
        
            FigDir = fullfile(pwd,'analysis','figure','reward_modulation_takeoff', ...
                sprintf('ephys%s',bname_ephys), ...
                trjName, sprintf('cluster%d',cluster_list(clu)));
            if isempty(dir(FigDir)); mkdir(FigDir); end
        
            % oth_dist_sorted = sort(oth_dist,'ascend');
            % social_flight = oth_dist_sorted <= 0.5; % social flight
            % nonsocial_flight = oth_dist_sorted > 0.5; % non-social flight
        
            % if sum(social_flight) > 1 && sum(nonsocial_flight) > 1
            if length(t_onset_sec_r1) > 3 && length(t_onset_sec_r2) > 3
        
                for unit = 1:n_unit
                    unit_id = unit_ids(unit);
                    depth = spikeData.good_units.depth(unit); % depth of unit
                    layers = Layers(strcmp({Layers.trjName},trjName) & strcmp({Layers.batid},bname_ephys)); % specify the layer
                    region = 'unknown';
                    for i = 1:length(layers.layer)
                        if depth > layers.layer(i).depth(1) && depth < layers.layer(i).depth(2)
                            region = layers.layer(i).name; % brain region of the current unit
                        end
                    end
            
                    sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6; % spike times of the current unit in sec
                   
                    % sptimes = alignSpikeTimes_OneUnit(sptimes_unit_sec,t_onset_sec,t_end_sec,t_mergin_sec); % extract spiketimes
                       
                    fig = figure('Visible','off');
                    % fig = figure;
                    set(fig, 'units','normalized','Position',[0.2 0.2 0.55 0.35]);
                    tl = tiledlayout(1,2,"TileSpacing",'compact');
                    title(tl, sprintf('%s, %s, unit %d, %s',bname_ephys,mtData.recName,unit_id,region),'interpreter','none')
                    
                    %%% FIGURE 1: raster plot around landings
                    nexttile
                     [sptimes_sec_r1,sptimes_sec_r2] = sortedRasterPlot_v3(sptimes_unit_sec, ...
                        t_onset_sec_r1,t_end_sec_r1,t_onset_sec_r2,t_end_sec_r2,t_margin);
                    hold on
                    line([0 0], [0 length(sptimes_sec_r1)+length(sptimes_sec_r2)],'Color',[150 150 150 255]/255,'LineWidth',2)
                    yregion(0,n_trj,'FaceColor',[150 150 150]/255)
                    yregion(0,sum(idx_r1),'FaceColor',[200 100 100]/255)
                    hold off
                                
                    xlim(t_margin)
                                
                    xlabel('Time from takeoff (sec)')
                    ylabel('Flights')
                    set(gca,'FontSize',14)
                    
                    %%% FIGURE 2: firing rate
                    nexttile
            
                    sp_bin = 0.1; % bin size (sec)
                    filt_smp = 5; % number of samples to smooth
                    sp_edges = t_margin(1)-sp_bin/2:sp_bin:t_margin(2)+sp_bin/2; % bin sptimes into 0.1ms bins

                    [fr_r1,sem_r1] = compute_psth(sptimes_sec_r1,sp_edges,sp_bin,filt_smp);
                    [fr_r2,sem_r2] = compute_psth(sptimes_sec_r2,sp_edges,sp_bin,filt_smp);
                    % N_r1 = zeros(length(sptimes_sec_r1),length(sp_edges)-1);
                    % N_r2 = zeros(length(sptimes_sec_r2),length(sp_edges)-1);
                    % 
                    % for ff = 1:length(sptimes_sec_r1)
                    %     [N,~] = histcounts(sptimes_sec_r1(ff).sptimes,sp_edges);
                    %     N_r1(ff,:) = N;
                    % end
                    % N_r1 = smoothdata(N_r1,2,'gaussian',filt_smp);
                    % sem_r1 = std(N_r1,0,1) / sqrt(size(N_r1,1));
                    % fr_r1 = mean(N_r1,1);
                    % % fr_social = smoothdata(fr_social,'gaussian',filt_smp);
                    % for ff = 1:length(sptimes_sec_r2)
                    %     [N,~] = histcounts(sptimes_sec_r2(ff).sptimes,sp_edges);
                    %     N_r2(ff,:) = N;
                    % end
                    % N_r2 = smoothdata(N_r2,2,'gaussian',filt_smp);
                    % sem_r2 = std(N_r2,0,1) / sqrt(size(N_r2,1));
                    % fr_r2 = mean(N_r2,1);
                    % fr_nonsocial = smoothdata(fr_nonsocial,'gaussian',filt_smp);
                    y_max = max(max(fr_r1+sem_r1),max(fr_r2 + sem_r2)) * 1.5; 
                                
            
                    % figure
                    hold on
                    % social curve
                    x = t_margin(1):sp_bin:t_margin(2);
                    % plot(x,fr_social,'r')
                    boundedline(x,fr_r1,sem_r1,'r','alpha','transparency', 0.1)
                               
                    % non-social curve
                    % plot(x,fr_nonsocial,'k')
                    boundedline(x,fr_r2,sem_r2,'k','alpha','transparency', 0.1)
                    line([0 0], [0 y_max],'Color',[150 150 150 255]/255,'LineWidth',2)
                    hold off
                    if y_max ~= 0
                        ylim([0 y_max])
                    end
                    xlabel('Time from takeoff (sec)')
            
                    set(gca,'FontSize',14)
                    
                    % Save
                    saveas(fig,fullfile(FigDir,sprintf('Unit%d.jpg',unit_id)))
                    close(fig)
                end
            end
        end
    end
end