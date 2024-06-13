function plot_social_modulation_land(bname_ephys,trjName,Flights,TrjCluster,spikeData,mtData,Layers)
% Specify the ephys bat
    bb = find(strcmp(bname_ephys,mtData.bname));
    
    % Load behavioral and neural features
    n_bats = length(mtData.bname); % number of bats
    oth_ids = [1:bb-1,bb+1:n_bats]; % index for other conspecifics
    t_land_sec = Flights(bb).t_land_sec;
    n_land = length(t_land_sec); % number of landings
    cluster_idx = TrjCluster(bb).clu_idx; % flight cluster index
    cluster_list = TrjCluster(bb).clu_list; % unique clusters
    cluster_count = TrjCluster(bb).clu_sz; % number of flights in each cluster
    n_cluster = length(cluster_count); % number of clusters
    n_unit = spikeData.numGoodUnits; % number of single units
    unit_ids = spikeData.good_units.cluster_id;
    
    
    t_margin = [-1 1]; % time window to extract spike times around landings (sec)
    
    
    for clu = 1:n_cluster
    % for clu = 29:29
        n_trj = cluster_count(clu); % number of trajectory for the current cluster
        
        if n_trj >= 3 && cluster_list(clu) ~= -1
            idx = cluster_idx == cluster_list(clu);
            t_onset_sec = t_land_sec(idx,1);
            t_end_sec = t_land_sec(idx,1);
            for k = 1:n_bats
                if k == 1 % plot for nearest neighbor bat
                    oth_dist = Flights(bb).closest_oth_dist(idx);  % distance to the nearest neighbor bat at landing
                    FigDir = fullfile(pwd,'analysis','figure','social_modulation_landing', ...
                        sprintf('ephys%s',bname_ephys), ...
                        trjName, sprintf('cluster%d',cluster_list(clu)), ...
                        'anyNN');
                    
                else
                    oth_id = oth_ids(k-1);
                    oth_dist = Flights(bb).dist_oth(idx,oth_id);
                    FigDir = fullfile(pwd,'analysis','figure','social_modulation_landing', ...
                        sprintf('ephys%s',bname_ephys), ...
                        trjName, sprintf('cluster%d',cluster_list(clu)), ...
                        sprintf('To%s',mtData.bname{oth_id}));

                end
                oth_dist_sorted = sort(oth_dist,'ascend');
                social_flight = oth_dist_sorted <= 0.5; % social flight
                nonsocial_flight = oth_dist_sorted > 0.5; % non-social flight
    
                if sum(social_flight) > 1 && sum(nonsocial_flight) > 1
                    if isempty(dir(FigDir)); mkdir(FigDir); end
        
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
                        set(fig, 'units','normalized','Position',[0.2 0.2 0.65 0.35]);
                        tl = tiledlayout(1,3,"TileSpacing",'compact');
                        title(tl, sprintf('%s, %s, unit %d, %s',bname_ephys,mtData.recName,unit_id,region),'interpreter','none')
            
                        %%% FIGURE 1: raster plot around landings
                        nexttile
                        sptimes_sec_sorted = sortedRasterPlot_v2(sptimes_unit_sec,t_onset_sec,t_end_sec,t_margin,oth_dist);
            
                        hold on
                        line([0 0], [0 length(sptimes_sec_sorted)],'Color',[150 150 150 255]/255,'LineWidth',2)
                        yregion(0,n_trj,'FaceColor',[150 150 150]/255)
                        yregion(0,sum(social_flight),'FaceColor',[200 100 100]/255)
                        hold off
                        
                        xlim(t_margin)
                        
                        xlabel('Time from landing (sec)')
                        ylabel('Flights')
                        set(gca,'FontSize',14)
            
                        %%% FIGURE 2: firing rate
                        nexttile
          
       
                        social_sptimes = sptimes_sec_sorted(social_flight);
                        nonsocial_sptimes = sptimes_sec_sorted(nonsocial_flight);
                        sp_bin = 0.1; % bin size (sec)
                        filt_smp = 5; % number of samples to smooth
                        sp_edges = t_margin(1)-sp_bin/2:sp_bin:t_margin(2)+sp_bin/2; % bin sptimes into 0.1ms bins

                        % Social flight
                        [fr_social,sem_social] = compute_psth(social_sptimes,sp_edges,sp_bin,filt_smp);
                        % 
                        % N_social = zeros(length(social_sptimes),length(sp_edges)-1);
                        % for ff = 1:length(social_sptimes)
                        %     [N,~] = histcounts(social_sptimes(ff).sptimes,sp_edges);
                        %     N_social(ff,:) = N;
                        % end
                        % N_social = smoothdata(N_social,2,'gaussian',filt_smp);
                        % sem_social = std(N_social,0,1) / sqrt(size(N_social,1));
                        % fr_social = mean(N_social,1);
                        % fr_social = smoothdata(fr_social,'gaussian',filt_smp);

                        % Non-social flight
                        [fr_nonsocial,sem_nonsocial] = compute_psth(nonsocial_sptimes,sp_edges,sp_bin,filt_smp);
                        % N_nonsocial = zeros(length(nonsocial_sptimes),length(sp_edges)-1);
                        % for ff = 1:length(nonsocial_sptimes)
                        %     [N,~] = histcounts(nonsocial_sptimes(ff).sptimes,sp_edges);
                        %     N_nonsocial(ff,:) = N;
                        % end
                        % N_nonsocial = smoothdata(N_nonsocial,2,'gaussian',filt_smp);
                        % sem_nonsocial = std(N_nonsocial,0,1) / sqrt(size(N_nonsocial,1));
                        % fr_nonsocial = mean(N_nonsocial,1);
                        % fr_nonsocial = smoothdata(fr_nonsocial,'gaussian',filt_smp);
                        y_max = max(max(fr_social+sem_social),max(fr_nonsocial + sem_nonsocial)) * 1.5; 
                        
    
                        % figure
                        hold on
                        % social curve
                        x = t_margin(1):sp_bin:t_margin(2);
                        % plot(x,fr_social,'r')
                        boundedline(x,fr_social,sem_social,'r','alpha','transparency', 0.1)
                       
                        % non-social curve
                        % plot(x,fr_nonsocial,'k')
                        boundedline(x,fr_nonsocial,sem_nonsocial,'k','alpha','transparency', 0.1)
                        line([0 0], [0 y_max],'Color',[150 150 150 255]/255,'LineWidth',2)
                        hold off
                        if y_max ~= 0
                            ylim([0 y_max])
                        end
                        xlabel('Time from landing (sec)')
    
                        set(gca,'FontSize',14)
            
                        %%% FIGURE 3: NNbat distance
                        nexttile
                        plot(oth_dist_sorted,1:length(oth_dist_sorted),'k')
    
                        hold on
                        yregion(0,n_trj,'FaceColor',[150 150 150]/255)
                        yregion(0,sum(social_flight),'FaceColor',[200 100 100]/255)
                        hold off
                        xlabel('Distance (m)')
                        ylabel('Flights')
    
                        set(gca,'FontSize',14)
    
                        % save
    
                        saveas(fig,fullfile(FigDir,sprintf('Unit%d.jpg',unit_id)))
                        close(fig)
                    end
                end
            end
    
        end     
    end