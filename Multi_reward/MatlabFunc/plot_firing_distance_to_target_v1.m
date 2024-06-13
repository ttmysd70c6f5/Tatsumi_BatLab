function plot_firing_distance_to_target_v1(bname_ephys,trjName,Flights,TrjCluster,spikeData,Layers,mtData)

% Specify the ephys bat
bb = find(strcmp(bname_ephys,mtData.bname));

% Load behavioral and neural features
n_bats = length(mtData.bname); % number of bats
oth_ids = [1:bb-1,bb+1:n_bats]; % index for other conspecifics
t_flight_sec = Flights(bb).t_flight_sec;
% n_land = length(t_land_sec); % number of landings
cluster_idx = TrjCluster(bb).clu_idx; % flight cluster
cluster_list = TrjCluster(bb).clu_list; % idx of clusters
cluster_sz = TrjCluster(bb).clu_sz; % number of flights in each cluster
n_cluster = length(cluster_sz); % number of clusters
n_unit = spikeData.numGoodUnits; % number of single units
unit_ids = spikeData.good_units.cluster_id;
    
t_margin = [0 0]; % time window to extract spike times around landings (sec)
    
    
for clu = 1:n_cluster
    n_trj = cluster_sz(clu); % number of trajectory for the current cluster    
    
    % if cluster_idx(clu) ~= -1
    if n_trj >= 3 && cluster_list(clu) ~= -1
        idx = cluster_idx == cluster_list(clu);
        t_onset_sec = t_flight_sec(idx,1);
        t_end_sec = t_flight_sec(idx,2);
        t_flights = Flights(bb).t_flight(idx);
        current_flights = TrjCluster(bb).r(idx);
        % for ll = 1:length(current_flights)
        %     current_flights{ll} = current_flights{ll}(:,:,bb);
        % end

        for k = 1:n_bats
            if k == 1 % plot for nearest neighbor bat
                oth_dist = Flights(bb).closest_oth_dist(idx);  % distance to the nearest neighbor bat at landing
                FigDir = fullfile(pwd,'analysis','figure','firing_distance_to_goal', ...
                        sprintf('ephys%s',bname_ephys), ...
                        trjName, sprintf('cluster%d',cluster_list(clu)), ...
                        'anyNN');
                    
            else
                oth_id = oth_ids(k-1);
                oth_dist = Flights(bb).dist_oth(idx,oth_id);
                FigDir = fullfile(pwd,'analysis','figure','firing_distance_to_goal', ...
                    sprintf('ephys%s',bname_ephys), ...
                    trjName, sprintf('cluster%d',cluster_list(clu)), ...
                    sprintf('To%s',mtData.bname{oth_id}));
            end

            oth_dist_sorted = sort(oth_dist,'ascend');
            social_flight = oth_dist_sorted <= 0.3; % social flight
            nonsocial_flight = oth_dist_sorted > 0.3; % non-social flight
    
            if sum(social_flight) > 3 && sum(nonsocial_flight) > 3
                if isempty(dir(FigDir)); mkdir(FigDir); end
    
                n_land = length(t_onset_sec);
    
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
    
                    % extract spike timings around the target behavioral features
                    sptimes = alignSpikeTimes_OneUnit(sptimes_unit_sec,t_onset_sec,t_end_sec,t_margin);
                    [splocs,dist_to_goal] = interpSplocs(t_flights,current_flights,sptimes);   
        
                    n_trials = length(t_onset_sec); % number of trials
                    % longest_trial = max(t_end_sec-t_onset_sec); % the longest trial length in sec
                    
                    % sort trials according to the event length
                    max_dist_goal = zeros(length(dist_to_goal),1);
                    for ll = 1:length(dist_to_goal)
                        max_dist_goal(ll) = max(dist_to_goal(ll).dist);
                    end
                    [~,idx_sort] = sort(oth_dist,'ascend'); % sort trials
                    max_dist_goal = max_dist_goal(idx_sort);
                    % t_end = t_end(idx_sort);
                    splocs_sorted = splocs(idx_sort);
    
                    %%% FIGURE
                    fig = figure('Visible','off');
                    set(fig, 'units','normalized','Position',[0.2 0.2 0.65 0.35]);
                    tl = tiledlayout(1,3,"TileSpacing",'compact');
                    title(tl, sprintf('%s, %s, unit %d, %s',bname_ephys,mtData.recName,unit_id,region),'interpreter','none')
       
                    %%% FIGURE 1: raster plot during the target behavior
                    % fig = figure('Visible','off');
                    nexttile
                    % figure
                    hold on
                    for ll = 1:n_trials
                        spX_trial = splocs_sorted(ll).splocs; % spike location in the current trial
                        % spY_trial = repmat(ll,size(spX_trial)); % location of spikes
                        for ss = 1:length(spX_trial)
                            line([spX_trial(ss) spX_trial(ss)], [ll-1 ll],'Color','k')
                            % plot(spX_trial,spY_trial,'.k','MarkerSize', 5);
                        end
                        line([0 0], [ll-1 ll],'Color',[200 100 100 255]/255,'LineWidth',1) % onset of trial
                        line([max_dist_goal(ll) max_dist_goal(ll)],[ll-1 ll],'Color',[200 100 100 255]/255,'LineWidth',1) % end of trial
                    end
                    yregion(0,n_trj,'FaceColor',[150 150 150]/255)
                    yregion(0,sum(social_flight),'FaceColor',[200 100 100]/255)
                    hold off

                    % title(sprintf('%s, %s, unit %d %s',bname_ephys,recDate,unit_id,region),'Interpreter','none')
                    xlabel('Distance to goal (m)')
                    ylabel('flight')
                    xlim([0, ceil(max(max_dist_goal))])
                    set(gca,'XDir','reverse','FontSize',14)


                    %%% FIGURE 2: PSTH
                    nexttile

                    [~,idx_sorted] = sort(oth_dist,'ascend'); % sort trials
                    splocs_sorted = splocs(idx_sorted);
                    social_splocs = splocs_sorted(social_flight);
                    nonsocial_splocs = splocs_sorted(nonsocial_flight);

                    sp_bin = 0.1; % bin size (m)
                    filt_smp = 10; % number of samples to smooth
                    sp_edges = 0:sp_bin:ceil(max(max_dist_goal));

                    % N_all = zeros(length(splocs_sorted),length(sp_edges)-1);
                    % for ll = 1:n_trials
                    %     [N_all(ll,:),~] = histcounts(splocs_sorted(ll).splocs,sp_edges);
                    % end
                    % N_all = smoothdata(N_all,2,'movmean',filt_smp) / sp_bin;
                    % sem_all = std(N_all,0,1)  / sqrt(size(N_all,1));
                    % fr_all = mean(N_all,1);
                    
                    N_social = zeros(length(social_splocs),length(sp_edges)-1);
                    for ll = 1:size(N_social,1)
                        [N_social(ll,:),~] = histcounts(social_splocs(ll).splocs,sp_edges);
                    end
                    N_social = smoothdata(N_social,2,'movmean',filt_smp);
                    sem_social = std(N_social,0,1)  / sqrt(size(N_social,1));
                    fr_social = mean(N_social,1);

                    N_nonsocial = zeros(length(nonsocial_splocs),length(sp_edges)-1);
                    for ll = 1:size(N_nonsocial,1)
                        [N_nonsocial(ll,:),~] = histcounts(nonsocial_splocs(ll).splocs,sp_edges);
                    end
                    N_nonsocial = smoothdata(N_nonsocial,2,'movmean',filt_smp);
                    sem_nonsocial = std(N_nonsocial,0,1)  / sqrt(size(N_nonsocial,1));
                    fr_nonsocial = mean(N_nonsocial,1);


                    % y_max = max([max(fr_all + sem_all), max(fr_social+sem_social),max(fr_nonsocial + sem_nonsocial)]) * 1.5;
                    y_max = max([max(fr_social+sem_social),max(fr_nonsocial + sem_nonsocial)]) * 1.5;


                    
                    ylabel('Firing rate')
                    xlabel('Distance to goal (m)')
                    xlim([0, ceil(max(max_dist_goal))])
                    set(gca,'XDir','reverse','FontSize',14)
    
                    % figure
                    t_fig = 0+sp_bin/2:sp_bin:ceil(max(max_dist_goal))-sp_bin/2;
                    hold on
                    % plot(t_fig,fr_all,'k')
                    boundedline(t_fig,fr_social,sem_social,'r','alpha','transparency', 0.1)
                    boundedline(t_fig,fr_nonsocial,sem_nonsocial,'k','alpha','transparency', 0.1)

                    % line([0 0], [0 y_max],'Color',[150 150 150 255]/255,'LineWidth',2)
                    hold off
                    if y_max ~= 0
                        ylim([0 y_max])
                    end
                    xlabel('Distance to goal (m)')

                    set(gca,'FontSize',14)
            

                    %%% FIGURE 3: NNbat distance
                    nexttile
                    plot(oth_dist_sorted,1:length(oth_dist_sorted),'k')

                    hold on
                    yregion(0,n_trj,'FaceColor',[150 150 150]/255)
                    yregion(0,sum(social_flight),'FaceColor',[200 100 100]/255)
                    hold off
                    xlabel('Distance to NN (m)')
                    ylabel('Flights')

                    set(gca,'FontSize',14)

                    % save
                    saveas(fig,fullfile(FigDir,sprintf('firing_distance_to_goal_unit%d.jpg',unit_id)))
                    close(fig)
                end
            end
        end

    end     
end

% % visualize firing for all flights for the same trajectory clusters
% for clu = 1:n_cluster
%     % if cluster_idx(clu) ~= -1
%     if cluster_list(clu) ~= -1
%         idx = cluster_idx == cluster_list(clu);
%         t_onset_sec = t_flight_sec(idx,1);
%         t_end_sec = t_flight_sec(idx,2);
%         t_flights = Flights(bb).t_flight(idx);
%         current_flights = TrjCluster(bb).r(idx);
% %         for ll = 1:length(current_flights)
% %             current_flights{ll} = current_flights{ll}(:,:,bb);
% %         end
% 
%         FigDir = fullfile(mtData.rootdir,'analysis','figure','firing_distance_to_goal', ...
%             sprintf('ephys%s',bname_ephys), ...
%             trjName, sprintf('cluster%d',cluster_list(clu)), 'all');
%         if isempty(dir(FigDir)); mkdir(FigDir); end
% 
%         n_land = length(t_onset_sec);
% 
%         for unit = 1:n_unit
%             unit_id = unit_ids(unit);
%             depth = spikeData.good_units.depth(unit); % depth of unit
%             layers = Layers(strcmp({Layers.trjName},trjName) & strcmp({Layers.batid},bname_ephys)); % specify the layer
%             region = 'unknown';
%             for i = 1:length(layers.layer)
%                 if depth > layers.layer(i).depth(1) && depth < layers.layer(i).depth(2)
%                     region = layers.layer(i).name; % brain region of the current unit
%                 end
%             end
% 
%             sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6; % spike times of the current unit in sec  
% 
%             % extract spike timings around the target behavioral features
%             sptimes = alignSpikeTimes_OneUnit(sptimes_unit_sec,t_onset_sec,t_end_sec,t_margin);
%             [splocs,dist_to_goal] = interpSplocs(t_flights,current_flights,sptimes);   
% 
%             n_trials = length(t_onset_sec); % number of trials
%             % longest_trial = max(t_end_sec-t_onset_sec); % the longest trial length in sec
% 
%             % sort trials according to the event length
%             max_dist_goal = zeros(length(dist_to_goal),1);
%             for ll = 1:length(dist_to_goal)
%                 max_dist_goal(ll) = max(dist_to_goal(ll).dist);
%             end
%             [~,idx_sort] = sort(max_dist_goal,'ascend'); % sort trials
%             max_dist_goal = max_dist_goal(idx_sort);
%             % t_end = t_end(idx_sort);
%             splocs_sorted = splocs(idx_sort);
% 
%             %%% FIGURE
%             fig = figure('Visible','off');
%             set(fig, 'units','normalized','Position',[0.2 0.2 0.55 0.35]);
%             tl = tiledlayout(1,2,"TileSpacing",'compact');
%             title(tl, sprintf('%s, %s, unit %d, %s',bname_ephys,recName,unit_id,region),'interpreter','none')
% 
%             %%% FIGURE 1: raster plot during the target behavior
%             % fig = figure('Visible','off');
%             nexttile
%             % figure
%             hold on
%             for ll = 1:n_trials
%                 spX_trial = splocs_sorted(ll).splocs; % spike location in the current trial
%                 % spY_trial = repmat(ll,size(spX_trial)); % location of spikes
%                 for ss = 1:length(spX_trial)
%                     line([spX_trial(ss) spX_trial(ss)], [ll-1 ll],'Color','k')
%                     % plot(spX_trial,spY_trial,'.k','MarkerSize', 5);
%                 end
%                 line([0 0], [ll-1 ll],'Color',[200 100 100 255]/255,'LineWidth',1) % onset of trial
%                 line([max_dist_goal(ll) max_dist_goal(ll)],[ll-1 ll],'Color',[200 100 100 255]/255,'LineWidth',1) % end of trial
%             end
%             % yregion(0,n_trj,'FaceColor',[150 150 150]/255)
%             % yregion(0,sum(social_flight),'FaceColor',[200 100 100]/255)
%             hold off
%             xlabel('Distance to goal (m)')
%             ylabel('flight')
%             xlim([0, ceil(max(max_dist_goal))])
%             set(gca,'XDir','reverse','FontSize',14)
% 
%             %%% FIGURE 2: PSTH
%             nexttile
% 
%             sp_bin = 0.1; % bin size (m)
%             filt_smp = 10; % number of samples to smooth
%             sp_edges = 0:sp_bin:ceil(max(max_dist_goal));
% 
%             N_all = zeros(length(splocs_sorted),length(sp_edges)-1);
%             for ll = 1:n_trials
%                 [N_all(ll,:),~] = histcounts(splocs_sorted(ll).splocs,sp_edges);
%             end
%             N_all = smoothdata(N_all,2,'gaussian',filt_smp);
%             sem_all = std(N_all,0,1) / sqrt(size(N_all,1));
%             fr_all = mean(N_all,1);
% 
%             y_max = max(fr_all + sem_all) * 1.5; 
% 
% 
% 
%             ylabel('Firing rate')
%             xlabel('Distance to goal (m)')
%             set(gca,'XDir','reverse','FontSize',14)
% 
%             % figure
%             t_fig = 0+sp_bin/2:sp_bin:ceil(max(max_dist_goal))-sp_bin/2;
%             hold on
%             plot(t_fig,fr_all,'k')
% 
%             % line([0 0], [0 y_max],'Color',[150 150 150 255]/255,'LineWidth',2)
%             hold off
%             if y_max ~= 0
%                 ylim([0 y_max])
%             end
%             xlabel('Distance to goal (m)')
%             xlim([0, ceil(max(max_dist_goal))])
% 
%             set(gca,'FontSize',14)
% 
% 
%             % save
%             saveas(fig,fullfile(FigDir,sprintf('firing_dist_to_goal_unit%d.jpg',unit_id)))
%             close(fig)
%         end
%     end
% end     