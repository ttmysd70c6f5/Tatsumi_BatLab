function genPlot_social_mod_oth_takeoff_v2(BhvData,Flights,TrjCluster,spikeData,mtData,FigDir)

% Check the distribution of self-accelerometer signals at others' landing
% Ephys bat
bb = find(strcmp(mtData.bname_ephys,mtData.bname));

% Behavioral variables
n_bat = BhvData.n_tags;
t = BhvData.t; % time
bnames = mtData.bname;
blanding = 1 - BhvData.bflying; % 1 if at resting sites
n_flight = Flights(bb).n_flight; % # of flights
t_land_idx = Flights(bb).t_land_idx; % time index of landings
d_oth_rest = Flights(bb).dist_oth_rest;
t_lands = Flights(bb).t_land;
ac_land = Flights(bb).a_land;

% Parameters
th_social = 0.5; % threshold for social flight (m)
t_margin = [-1.5 1.5]; % window range to evaluate social modualtion (sec)

for bbb = [1:bb-1,bb+1:n_bat]
    % bbb = 8;
    
    % % Identify social flights
    % n_flight = Flights(bbb).n_flight; % number of flights
    % t_oth_toff_sec = Flights(bbb).t_flight_sec(:,1); % Time of landings of the target bat
    % % d_oth_land = Flights(bbb).dist_oth(:,bb); % distance to the ephys bat at landings of the target bat
    % r_toff = cell2mat(cellfun(@(x) x(1,:,[bb,bbb]), Flights(bbb).r_flight, 'UniformOutput',false));
    % d_toff = vecnorm(diff(r_toff,1,3),2,2);
    % is_social = d_toff < th_social;

    t_oth_toff_sec = []; % conspecific land timings
    is_social = [];
    ac = [];
    for ll = 1:n_flight
        is_ctoff = diff([blanding(t_land_idx(ll,1),bbb); blanding(t_land_idx(ll,1):t_land_idx(ll,2),bbb)]) == -1;
        t_oth_toff_sec = [t_oth_toff_sec; t_lands{ll}(is_ctoff)];
        is_social = [is_social; d_oth_rest{ll}(is_ctoff,bbb) < th_social];
        ac = [ac; ac_land{ll}(is_ctoff,bb)];
    end
    is_social = logical(is_social);
    
    % Exclude the flight at very biggining or ending
    is_include = t_oth_toff_sec > t(1)-t_margin(1) & t_oth_toff_sec < t(end)-t_margin(2);
    
    % Exclude the flight without ephys recording
    t_ephys_start = min(spikeData.allGlobalSpikeTimes_usec/1e6);
    t_ephys_end = max(spikeData.allGlobalSpikeTimes_usec/1e6);
    is_include2 = t_oth_toff_sec > t_ephys_start & t_oth_toff_sec < t_ephys_end;
    
    %% Distribution of acceleration for social/non-social flights
    % ac = cell2mat(cellfun(@(x) x(1,bb), Flights(bbb).a_flight, 'UniformOutput',false)); % acceleration (x) of self at takeoff
    n_trjclu = TrjCluster(bbb).n_clu;
    
    %%% FIGURE
    fig = figure('Visible','off');
    % fig = figure;
    set(fig, 'units','normalized','Position',[0 0 1 1]);
    tl = tiledlayout(3,6,"TileSpacing",'compact');
    title(tl,sprintf('Bat %s',mtData.bname{bbb}))
    for clu = 1:n_trjclu
        nexttile
        clu_now = TrjCluster(bbb).clu_idx == TrjCluster(bbb).clu_list(clu); % logical vector for the current cluster
        hold on
        % social
        idx_social = is_social&is_include&is_include2&clu_now;
        histogram(ac(idx_social),floor(min(ac)):0.25:ceil(max(ac)),'FaceColor','r','FaceAlpha',0.3)
        % non-social
        idx_asocial = ~is_social&is_include&is_include2&clu_now;
        histogram(ac(idx_asocial),floor(min(ac)):0.25:ceil(max(ac)),'FaceColor','k','FaceAlpha',0.3)
        hold off
    
        title(sprintf('Clu%d: Social:%d,Non-social:%d',TrjCluster(bbb).clu_list(clu),sum(idx_social),sum(idx_asocial)));
        xlabel('')
        set(gca,'FontSize',14)
    end
    
    saveDir = fullfile(FigDir,'ac_adjust');
    if isempty(dir(saveDir)); mkdir(saveDir); end
    saveas(fig,fullfile(saveDir,sprintf('ac_before_adjust_%s.png',bnames{bbb})))
    
    %% Adjust the acceleration distribution
    seed = 1234;
    
    %%% FIGURE
    fig = figure('Visible','off');
    % fig = figure;
    set(fig, 'units','normalized','Position',[0 0 1 1]);
    tl = tiledlayout(3,6,"TileSpacing",'compact');
    title(tl,sprintf('Bat %s',mtData.bname{bbb}))
    
    % Array to save the index of downsampled flights
    social_adjusted = NaN(1,n_flight); asocial_adjusted = NaN(1,n_flight);
    for clu = 1:n_trjclu   
        % Identify social/asocial flights
        clu_now = TrjCluster(bbb).clu_idx == TrjCluster(bbb).clu_list(clu); % logical vector for the current cluster
        idx_social = find(is_social&is_include&is_include2&clu_now);
        idx_asocial = find(~is_social&is_include&is_include2&clu_now);
    
        % Adjust distribution
        bin = 0.25;
        [ac_social_adjusted,ac_asocial_adjusted,social_adjusted_idx,asocial_adjusted_idx] = adjust_distribution(ac(idx_social),ac(idx_asocial),bin,seed);
        % Flight index after adjustment
        social_adjusted(idx_social(social_adjusted_idx)) = TrjCluster(bbb).clu_list(clu);
        asocial_adjusted(idx_asocial(asocial_adjusted_idx)) = TrjCluster(bbb).clu_list(clu);
    
        % Plot adjsuted distributions
        nexttile
        
        hold on
        histogram(ac_social_adjusted,floor(min(ac)):0.25:ceil(max(ac)),'FaceColor','r','FaceAlpha',0.3)    
        histogram(ac_asocial_adjusted,floor(min(ac)):0.25:ceil(max(ac)),'FaceColor','k','FaceAlpha',0.3)
        hold off
    
        title(sprintf('Clu%d: Social:%d,Non-social:%d',TrjCluster(bbb).clu_list(clu),length(ac_social_adjusted),length(ac_asocial_adjusted)));
        xlabel('')
        set(gca,'FontSize',14)
    end
    
    saveDir = fullfile(FigDir,'ac_adjust');
    if isempty(dir(saveDir)); mkdir(saveDir); end
    saveas(fig,fullfile(saveDir,sprintf('ac_after_adjust_%s.png',bnames{bbb})))
    
    %% Plot firing curve and raster plots for each cluster
    % Parameters
    w_step = 0.1; % step size of window for optimizing time bin (sec)
    w_sz = 0.5; % window size for optimizing time bin (sec)
    
    t_bin = t_margin(1):w_step:t_margin(2);
    pts_margin = 5; % marginal window to compute moving mean/sum
    
    % Load behavioral features
    % r = BhvData.r(:,:,bb); % position
    % v = BhvData.v_abs(:,bb); % velocity
    a_abs = BhvData.ac_abs(:,bb); % acceleration
    t = BhvData.t; % time
    % n_bats = BhvData.n_tags; % # of bats
    % blanding = 1 - BhvData.bflying; % 1 if at resting sites
    fs_bhv = BhvData.fs; % sampling frequency for behavioral data
    % n_flight_clu = Flights(bb).n_flight; % # of flights
    % t_land_idx = Flights(bb).t_land_idx; % time index of landings
    % d_oth_rest = Flights(bb).dist_oth_rest;
    % t_lands = Flights(bb).t_land;
    

    % Load neural features
    % target_units =  spikeData.good_units.fr >= 0.2 & ~strcmp(spikeData.good_units.region,'unknown') & ~strcmp(spikeData.good_units.region,'VC'); % firing rate criteria
    target_units =  spikeData.good_units.fr >= 0.2; % firing rate criteria
    unit_ids = spikeData.good_units.cluster_id(target_units);
    sptimes_all_usec = spikeData.good_units.spikeTimes_usec(target_units);
    unit_depths = spikeData.good_units.depth(target_units);
    n_unit = sum(target_units);
    frs = spikeData.good_units.fr(target_units);
    regions = spikeData.good_units.region(target_units);

    mod_other_toff(bb).clu_id = TrjCluster(bbb).clu_list;
    mod_other_toff(bb).id = unit_ids;
    mod_other_toff(bb).depth = unit_depths;
    mod_other_toff(bb).region = regions;
    mod_other_toff(bb).fr = frs;
    
    gen_fig = true(1,n_trjclu);
    for clu = 1:n_trjclu 
        clu_id = TrjCluster(bbb).clu_list(clu);
    

        % 1. Check exclusion criteria
        n_flight_clu = sum(social_adjusted == clu_id);
        if n_flight_clu < 5 % minimum number of flights
            gen_fig(clu) = false; % false if the # of flights is insufficient
        else
            % social/asocial flights after adjustment
            t_social = t_oth_toff_sec(social_adjusted == clu_id);
            t_asocial = t_oth_toff_sec(asocial_adjusted == clu_id);
    
            % number of flights
            n1 = n_flight_clu; n2 = n_flight_clu; % number of flights for each category (after adjustment)
            n_event = n1 + n2; % total length of flights
    
            %%%%% Run for all units
            % Arrays to save outputs
            t_opt_bin = NaN(1,n_unit);
            p_ctest = NaN(1,n_unit);
            p_glm = NaN(1,n_unit);
    
            for unit = 1:n_unit
            % for unit = 34:34 % run for all units
                disp(sprintf('Unit%d...',unit))
                % Unit features
                unit_id = unit_ids(unit);
                fr = frs(unit);
                depth = unit_depths(unit); % depth of unit
                region = regions{unit}; % brain region
                sptimes_unit_sec = sptimes_all_usec{unit}/1e6; % spike times of the current unit in sec
    
                % 2. Determine the optimal time bin
                % Count spikes for each time bin
                spmat_social = NaN(n1,diff(t_margin)/w_step+1+pts_margin*2); % spike count for each time bin (# of flights, # of time bins)
                spmat_asocial = NaN(n2,diff(t_margin)/w_step+1+pts_margin*2); % spike count for each time bin (# of flights, # of time bins)
                t_edges_social = [];
                t_edges_asocial = [];
                for ff = 1:n1
                    t_edges_social(ff,:) = t_social(ff)+t_margin(1)-w_step/2-pts_margin*w_step:w_step:t_social(ff)+t_margin(2)+w_step/2+pts_margin*w_step;
                    spmat_social(ff,:) = histcounts(sptimes_unit_sec,t_edges_social(ff,:));
                end
                for ff = 1:n2
                    t_edges_asocial(ff,:) = t_asocial(ff)+t_margin(1)-w_step/2-pts_margin*w_step:w_step:t_asocial(ff)+t_margin(2)+w_step/2+pts_margin*w_step;
                    spmat_asocial(ff,:) = histcounts(sptimes_unit_sec,t_edges_asocial(ff,:));
                end
                % Poisson mean difference test
                stat_test_poisson = NaN(3,size(spmat_social,2)); % row 1 = # of spikes for social flights, row 2 = # of spikes for non-social flights, row 3 = p-values for e-test
                stat_test_poisson(1,:) = movsum(sum(spmat_social,1),w_sz/w_step);
                stat_test_poisson(2,:) = movsum(sum(spmat_asocial,1),w_sz/w_step);
                for i = 1:size(stat_test_poisson,2)
                    % stat_test_poisson(3,i) = pois_mean_e_test(e_test_poisson(1,i),e_test_poisson(2,i),'n1',n1,'n2',n2);
                    stat_test_poisson(3,i) = poisson_c_test(stat_test_poisson(1,i),stat_test_poisson(2,i),n1,n2);  
                end
                stat_test_poisson(3,:) = stat_test_poisson(3,:) * size(stat_test_poisson,2); % adjusted p-values by Bon-Ferroni correction
                % Compute firing rate and standard deviation for each bin
                fr1 = stat_test_poisson(1,:)/w_sz/n1; % firing rates for social flights
                fr2 = stat_test_poisson(2,:)/w_sz/n2; % firing rates for non-social flights
                sem1 = std(movsum(spmat_social,w_sz/w_step,2)/w_sz,0,1) / n1; % standard error of mean for sllinding window (social flights)
                sem2 = std(movsum(spmat_asocial,w_sz/w_step,2)/w_sz,0,1) / n2; % standard error of mean for sliding window (non-social flights)
    
                % Screening the optimal time bin with a 500-ms time window
                [~,opt_bin] = min(stat_test_poisson(3,pts_margin+1:end-pts_margin));
                t_opt_bin(unit) = t_bin(opt_bin);
                p_ctest(unit) = stat_test_poisson(3,opt_bin+pts_margin); % p-value of the optimal bin
    
                % 3. Response variable: # of spikes for optimal bin (500ms window)
                % Compute the total # of spikes for the optimal bin
                spcnt_social = movsum(spmat_social,w_sz/w_step,2); spcnt_social = spcnt_social(:,pts_margin+opt_bin);
                spcnt_asocial = movsum(spmat_asocial,w_sz/w_step,2); spcnt_asocial = spcnt_asocial(:,pts_margin+opt_bin);
                
                % 5. Stepwise GLM
                % Extract explanatory variables
                % idx_opt_bin = ismember(round(t,2),round(t_opt_bin,2));
                % [~,idx_opt_bin] = ismember(round(t_opt_bin(unit),2),round(t,2));
                % r_x = r(idx_opt_bin,1);
                % r_y = r(idx_opt_bin,2);
                % v_abs = v(idx_opt_bin);
                % is_fd = double(is_fd);
                ac_glm = [ac(social_adjusted == clu_id); ac(asocial_adjusted == clu_id)]; % accelerometer signals
                social_glm = [ones(n1,1);zeros(n2,1)];
                spcnt_glm = [spcnt_social;spcnt_asocial]+1e-10; % avoid log(0)
            
                T_vars = table(ac_glm,social_glm,spcnt_glm);
                % T_vars = table(r_x,r_y,v_abs,is_social,sp_counts);
                % T_vars = table(is_social,v_abs,r_x,sp_counts);
                % mdl = stepwiseglm(T_vars,'constant','Lower','linear','Distribution','poisson','Link','log')
                mdl = stepwiseglm(T_vars,'constant','Distribution','poisson','Link','log');
                coef_names = mdl.CoefficientNames(2:end);
                coef_pvalues = mdl.Coefficients.pValue(2:end);
                coef_sig = coef_names(coef_pvalues < 0.001); % significant coefficients
                coef_sig = join(coef_sig,',');
                if any(strcmp(coef_names,'social_glm'))
                    p_glm(unit) = coef_pvalues(strcmp(coef_names,'social_glm'));
                end
    
    
                % Remove the margin bins at both edges
                fr1 = fr1(:,pts_margin+1:end-pts_margin);
                fr2 = fr2(:,pts_margin+1:end-pts_margin);
                sem1 = sem1(pts_margin+1:end-pts_margin);
                sem2 = sem2(pts_margin+1:end-pts_margin);
                spmat_social = spmat_social(:,pts_margin+1:end-pts_margin);
                spmat_asocial = spmat_asocial(:,pts_margin+1:end-pts_margin);
                t_edges_social = t_edges_social(:,pts_margin+1:end-pts_margin);
                t_edges_asocial = t_edges_asocial(:,pts_margin+1:end-pts_margin);
    
                % Spike arrays for raster plots
                spmat2_social = NaN(n1,diff(t_margin)*fs_bhv+1);
                spmat2_asocial = NaN(n2,diff(t_margin)*fs_bhv+1);
                for ff = 1:n1
                    t_edges_2 = t_social(ff)+t_margin(1)-1/fs_bhv/2:1/fs_bhv:t_social(ff)+t_margin(2)+1/fs_bhv/2;
                    spmat2_social(ff,:) = histcounts(sptimes_unit_sec,t_edges_2);
                end
                for ff = 1:n2
                    t_edges_2 = t_asocial(ff)+t_margin(1)-1/fs_bhv/2:1/fs_bhv:t_asocial(ff)+t_margin(2)+1/fs_bhv/2;
                    spmat2_asocial(ff,:) = histcounts(sptimes_unit_sec,t_edges_2);
                end
                
                % Peri-event firing rate (social vs non-social)
                pefr_social = sum(spmat2_social,'all')/(n1*diff(t_margin));
                pefr_asocial = sum(spmat2_asocial,'all')/(n2*diff(t_margin));
    
    
                %%% FIGURE: Firing rate curve and raster plot
                %%% Initialize figure bundle
                fig = figure('visible','off');
                % fig = figure;
                set(fig,'units','normalize','Position',[0.2 0.2 0.3 0.7])
                tl = tiledlayout(2,1,'TileSpacing','tight');
                title_str = sprintf('Unit%d, %s, %dum, %.1fHz (Social: %.1fHz, Non-social: %.1fHz)',unit_id,region,depth,fr,pefr_social,pefr_asocial);
                title(tl,title_str,'Interpreter','None','FontSize',14)
                % Plot 1: Firing rate curve
                nexttile
                t_fig = t_margin(1):w_step:t_margin(2); % time vector for figure
               
                hold on
                % plot(t_fig,e_test_poisson(1,:)/w_sz/n1,'r','DisplayName','Social')
                % plot(t_fig,e_test_poisson(2,:)/w_sz/n2,'k','DisplayName','Non-social')
                boundedline(t_fig,fr1,sem1,'color',[250 50 50]/255,'alpha','transparency', 0.1)
                boundedline(t_fig,fr2,sem2,'color',[100 100 100]/255,'alpha','transparency', 0.1)
                hold off
                ylabel('Firing rate (Hz)')
                % xlabel('Time from takeoffs (sec)')
                xlim(t_margin)
                title_str = sprintf('Others'' takeoffs: p=%.3f \nGLM:%s',min(stat_test_poisson(3,:)),coef_sig{:});
                % title_str = sprintf('Others'' landings: \nGLM:%s',coef_sig{:});
                title(title_str,'Interpreter','None')
                set(gca,'FontSize',14,'XColor','none','TickDir','none')
            
            
                % Plot 2: raster plot
                nexttile
                t_bins = t_margin(1):1/fs_bhv:t_margin(2);
                A = [spmat2_social; spmat2_asocial]; % array for raster plot
                A(A==0) = nan; A(A>0) = 1; % convert to binary
                h = imagesc([t_bin(1) t_bin(end)],[0.5 n1+n2-0.5],A);
                set(gca,'YDir','normal','FontSize',16)
                set(h,'AlphaData',~isnan(A))
                colormap gray
    
                hold on
                yregion(0,n1,'FaceColor',[250 50 50]/255,'FaceAlpha',0.2)
                yregion(n1,n1+n2,'FaceColor',[100 100 100]/255,'FaceAlpha',0.2)
                hold off
                xlim(t_margin)
                ylim([0,n1+n2])
                xlabel('Time from others'' takeoffs (sec)')
                yticks([0.5 n1+n2-0.5]);
                yticklabels({1,n1+n2})
                set(gca,'TickDir','none','FontSize',14)
    
            
                % the arrows
                xO = -0.2;  
                yO = 0.5;
                x_fd1 = t_margin(1)-0.1;
                x_fd2 = x_fd1;
                y_fd1 = round(n1/2);
                y_fd2 = round(n1+n2/2);
                % patch(...
                %     [x_fd1+xO x_fd2+xO; x_fd1+xO x_fd1+xO; x_fd1 x_fd2], ...
                %     [y_fd1-yO y_fd2-yO; y_fd1+yO y_fd2+yO; y_fd1 y_fd2], 'k', 'clipping', 'off')
            
                % Texts
                text(x_fd1-0.1, y_fd1, 'Social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
                text(x_fd2-0.1, y_fd2, 'Non-social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
            
                %%% Save figure
                saveDir = fullfile(FigDir,sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName),region,bnames{bbb},sprintf('Clu%d',clu_id));
                if isempty(dir(saveDir)); mkdir(saveDir); end
                saveas(fig,fullfile(saveDir,sprintf('Unit%d.jpg',unit_id)))
                close(fig)
            end
            mod_other_toff(bb).t_opt(clu).t_opt = t_opt_bin; % t with most significant modulation
            mod_other_toff(bb).p_ctest(clu).p_ctest = p_ctest;
            mod_other_toff(bb).p_glm(clu).p_glm = p_glm;
        end
    end
    mod_other_toff(bb).gen_fig = gen_fig;
end
saveFile = fullfile(FigDir,'social_modulation_other_toff.mat');
save(saveFile,'mod_other_toff');


% %%
% % Parameters
% t_margin = [-1.5 1.5]; % window range to evaluate social modualtion (sec)
% th_social = 0.5; % threshold for social flight (m)
% w_step = 0.1; % step size of window for optimizing time bin (sec)
% w_sz = 0.5; % window size for optimizing time bin (sec)
% 
% t_bin = t_margin(1):w_step:t_margin(2);
% pts_margin = 5; % marginal window to compute moving mean/sum
% 
% % Load behavioral features
% r = BhvData.r(:,:,bb); % position
% v = BhvData.v_abs(:,bb); % velocity
% t = BhvData.t; % time
% n_bats = BhvData.n_tags; % # of bats
% blanding = 1 - BhvData.bflying; % 1 if at resting sites
% fs_bhv = BhvData.fs; % sampling frequency for behavioral data
% n_flight_clu = Flights(bb).n_flight; % # of flights
% t_land_idx = Flights(bb).t_land_idx; % time index of landings
% d_oth_rest = Flights(bb).dist_oth_rest;
% t_lands = Flights(bb).t_land;
% % is_social = Flights(bb).closest_oth_dist < th_social; % true if social flight
% 
% 
% % Extract the timing of conspecifics landing
% t_oth_land_sec = []; % conspecific land timings
% is_social = [];
% for bbb = [1:bb-1,bb+1:n_bats]
%     for ff = 1:n_flight_clu
%         idx_land = diff([blanding(t_land_idx(ff,1),bbb); blanding(t_land_idx(ff,1):t_land_idx(ff,2),bbb)]) == 1;
%         t_oth_land_sec = [t_oth_land_sec; t_lands{ff}(idx_land)];
%         is_social = [is_social; d_oth_rest{ff}(idx_land,bbb) < th_social];
%     end
% end
% is_social = logical(is_social);
% [t_oth_land_sec,sort_idx] = sort(t_oth_land_sec);
% is_social = is_social(sort_idx);
% 
% % Exclude the flight at very biggining or ending
% is_include = t_oth_land_sec > t(1)-t_margin(1) & t_oth_land_sec < t(end)-t_margin(2);
% t_oth_land_sec = t_oth_land_sec(is_include);
% is_social = is_social(is_include);
% 
% t_oth_land_sec = t_oth_land_sec; % time of target events
% n_event = length(t_oth_land_sec); % number of target events
% n1 = sum(is_social); n2 = sum(~is_social);
% 
% % 1. Check exclusion criteria
% if n_flight_clu < 20 % minimum number of flights
%     error('The number of flights is insufficient.')
% end
% % if n1 < 10 || n2 < 10 % minimum number of flights per category
% %     error('The number of flights per category is insufficient.')
% % end
% % target_units =  spikeData.good_units.fr >= 0.2; % firing rate criteria
% % target_units =  spikeData.good_units.fr >= 0.2 & strcmp(spikeData.good_units.region,'MEC'); % firing rate criteria
% % target_units =  spikeData.good_units.fr >= 0.2; % firing rate criteria
% target_units =  spikeData.good_units.fr >= 0.2 & ~strcmp(spikeData.good_units.region,'unknown') & ~strcmp(spikeData.good_units.region,'VC'); % firing rate criteria
% 
% % Load neural features
% unit_ids = spikeData.good_units.cluster_id(target_units);
% sptimes_all_usec = spikeData.good_units.spikeTimes_usec(target_units);
% unit_depths = spikeData.good_units.depth(target_units);
% n_unit = sum(target_units);
% frs = spikeData.good_units.fr(target_units);
% 
% % Exclude the behavioral data after the ephys recording stop
% t_ephys_start = min(spikeData.allGlobalSpikeTimes_usec/1e6);
% t_ephys_end = max(spikeData.allGlobalSpikeTimes_usec/1e6);
% is_include = t_oth_land_sec > t_ephys_start & t_oth_land_sec < t_ephys_end;
% t_oth_land_sec = t_oth_land_sec(is_include);
% is_social = is_social(is_include);
% n_event = length(t_oth_land_sec);
% n1 = sum(is_social);
% n2 = sum(~is_social);
% 
% clear mod_other_land
% 
% %%%%% Run for all units
% % for unit = 34:34 % run for all units
% for unit = 1:n_unit
% % for unit = [1:17,19:22,24:n_unit]
% % for unit =  34:n_unit
%     disp(sprintf('Unit%d',unit))
%     % Unit features
%     unit_id = unit_ids(unit);
%     fr = frs(unit);
%     depth = unit_depths(unit); % depth of unit
%     layers = Layers(strcmp({Layers.trjName},trjName) & strcmp({Layers.batid},bname_ephys)); % specify the layer
%     region = get_probe_region(layers,depth);
%     sptimes_unit_sec = sptimes_all_usec{unit}/1e6; % spike times of the current unit in sec
% 
%     mod_other_land(unit).id = unit_id;
%     mod_other_land(unit).depth = depth;
%     mod_other_land(unit).region = region;
%     mod_other_land(unit).fr = fr;
% 
%     % 2. Determine the optimal time bin
%     % Count spikes for each time bin
%     spmat_social = zeros(n_event,diff(t_margin)/w_step+1+pts_margin*2); % spike count for each time bin (# of flights, # of time bins)
%     spmat2_social = zeros(n_event,diff(t_margin)*fs_bhv+1);
%     t_edges_social = [];
%     for ff = 1:n_event
%         t_edges_social(ff,:) = t_oth_land_sec(ff)+t_margin(1)-w_step/2-pts_margin*w_step:w_step:t_oth_land_sec(ff)+t_margin(2)+w_step/2+pts_margin*w_step;
%         t_edges_2 = t_oth_land_sec(ff)+t_margin(1)-1/fs_bhv/2:1/fs_bhv:t_oth_land_sec(ff)+t_margin(2)+1/fs_bhv/2;
%         spmat_social(ff,:) = histcounts(sptimes_unit_sec,t_edges_social(ff,:));
%         spmat2_social(ff,:) = histcounts(sptimes_unit_sec,t_edges_2);
%     end
% 
%     % Screening with 500-ms time window
%     stat_test_poisson = NaN(3,size(spmat_social,2)); % row 1 = # of spikes for social flights, row 2 = # of spikes for non-social flights, row 3 = p-values for e-test
%     stat_test_poisson(1,:) = movsum(sum(spmat_social(is_social,:),1),w_sz/w_step);
%     stat_test_poisson(2,:) = movsum(sum(spmat_social(~is_social,:),1),w_sz/w_step);
%     for i = 1:size(stat_test_poisson,2)
%         % e_test_poisson(3,i) = pois_mean_diff_test(e_test_poisson(1,i),e_test_poisson(2,i),'n1',n1,'n2',n2);
%         stat_test_poisson(3,i) = poisson_c_test(stat_test_poisson(1,i),stat_test_poisson(2,i),n1,n2);  
%     end
%     stat_test_poisson(3,:) = stat_test_poisson(3,:) * size(stat_test_poisson,2); % Bon-Ferroni correction
% 
%     fr1 = stat_test_poisson(1,:)/w_sz/n1; % firing rates for social flights
%     fr2 = stat_test_poisson(2,:)/w_sz/n2; % firing rates for non-social flights
%     sem1 = std(movsum(spmat_social(is_social,:),w_sz/w_step,2)/w_sz,0,1) / n1; % standard error of mean for sllinding window (social flights)
%     sem2 = std(movsum(spmat_social(~is_social,:),w_sz/w_step,2)/w_sz,0,1) / n2; % standard error of mean for sliding window (non-social flights)
% 
%     fr1 = fr1(:,pts_margin+1:end-pts_margin);
%     fr2 = fr2(:,pts_margin+1:end-pts_margin);
%     sem1 = sem1(pts_margin+1:end-pts_margin);
%     sem2 = sem2(pts_margin+1:end-pts_margin);
%     spmat_social = spmat_social(:,pts_margin+1:end-pts_margin);
%     t_edges_social = t_edges_social(:,pts_margin+1:end-pts_margin);
% 
%     % optimal time bin
%     [~,opt_bin] = min(stat_test_poisson(3,pts_margin+1:end-pts_margin));
%     % opt_bin = 13;
%     t_opt_bin = t_bin(opt_bin);
%     mod_other_land(unit).t_opt = t_opt_bin; % t with most significant modulation
%     mod_other_land(unit).p_stat = stat_test_poisson(3,opt_bin+pts_margin);
% 
%     % Peri-event firing rate (social vs non-social)
%     pefr_social = sum(spmat2_social(is_social,:),'all')/(n1*diff(t_margin));
%     pefr_asocial = sum(spmat2_social(~is_social,:),'all')/(n2*diff(t_margin));
% 
%     mod_other_land(unit).pefr_social = pefr_social;
%     mod_other_land(unit).pefr_nonsocial = pefr_asocial;
% 
%     % 3. Response variable: # of spikes for optimal bin (500ms window)
%     % sp_counts = sum(sp_mat(:,opt_bin-2:opt_bin+2),2);
%     spcnt_social = movsum(spmat_social,w_sz/w_step,2); spcnt_social = spcnt_social(:,opt_bin);
% 
%     % 4. Extract explanatory variables
%     % idx_opt_bin = ismember(round(t,2),round(t_opt_bin,2));
%     [~,idx_opt_bin] = ismember(round(t_opt_bin,2),round(t,2));
%     % r_x = r(idx_opt_bin,1);
%     % r_y = r(idx_opt_bin,2);
%     % v_abs = v(idx_opt_bin);
%     % is_fd = double(is_fd);
% 
%     % 5. Stepwise GLM
%     T_vars = table(is_social,spcnt_social);
%     % T_vars = table(r_x,r_y,v_abs,is_social,sp_counts);
%     % T_vars = table(is_social,v_abs,r_x,sp_counts);
%     % mdl = stepwiseglm(T_vars,'constant','Lower','linear','Distribution','poisson','Link','log')
%     mdl = stepwiseglm(T_vars,'constant','Distribution','poisson','Link','log');
%     coef_names = mdl.CoefficientNames(2:end);
%     coef_pvalues = mdl.Coefficients.pValue(2:end);
%     coef_sig = coef_names(coef_pvalues < 0.01); % significant coefficients
%     coef_sig = join(coef_sig,',');
% 
%     mod_other_land(unit).p_glm = coef_pvalues(strcmp(coef_names,'is_social_1'));
% 
%     %%% FIGURE: Firing rate curve and raster plot
%     %%% Initialize figure bundle
%     fig = figure('visible','off');
%     % fig = figure;
%     set(fig,'units','normalize','Position',[0.2 0.2 0.3 0.7])
%     tl = tiledlayout(2,1,'TileSpacing','tight');
%     title_str = sprintf('Unit%d, %s, %dum, %.1fHz (Social: %.1fHz, Non-social: %.1fHz)',unit_id,region,depth,fr,pefr_social,pefr_asocial);
%     title(tl,title_str,'Interpreter','None','FontSize',14)
%     % Plot 1: Firing rate curve
%     nexttile
%     t_fig = t_margin(1):w_step:t_margin(2); % time vector for figure
% 
%     hold on
%     % plot(t_fig,e_test_poisson(1,:)/w_sz/n1,'r','DisplayName','Social')
%     % plot(t_fig,e_test_poisson(2,:)/w_sz/n2,'k','DisplayName','Non-social')
%     boundedline(t_fig,fr1,sem1,'color',[250 50 50]/255,'alpha','transparency', 0.1)
%     boundedline(t_fig,fr2,sem2,'color',[100 100 100]/255,'alpha','transparency', 0.1)
%     hold off
%     ylabel('Firing rate (Hz)')
%     % xlabel('Time from takeoffs (sec)')
%     xlim(t_margin)
%     title_str = sprintf('Others'' landings: p=%.3f \nGLM:%s',min(stat_test_poisson(3,:)),coef_sig{:});
%     % title_str = sprintf('Others'' landings: \nGLM:%s',coef_sig{:});
%     title(title_str,'Interpreter','None')
%     set(gca,'FontSize',14,'XColor','none','TickDir','none')
% 
% 
%     % Plot 2: raster plot
%     nexttile
%     t_bins = t_margin(1):1/fs_bhv:t_margin(2);
%     [~,sort_idx] = sort(is_social);
%     sp_mat_2_sort = spmat2_social(sort_idx,:);
%     t_fig = []; y_loc = [];
%     for ff = 1:size(sp_mat_2_sort,1)
%         x_loc = find(sp_mat_2_sort(ff,:));
%         t_fig = [t_fig, t_bins(x_loc)];
%         y_loc = [y_loc, repmat(ff,1,length(x_loc))];
%     end
%     t_fig = repmat(t_fig,2,1); y_loc = [y_loc; y_loc + 1];
%     hold on
%     line(t_fig,y_loc,'Color','k','LineWidth',1.3)
%     yregion(n2,n1+n2+1,'FaceColor',[250 50 50]/255,'FaceAlpha',0.2)
%     yregion(0,n2,'FaceColor',[100 100 100]/255,'FaceAlpha',0.2)
%     hold off
%     xlim(t_margin)
%     ylim([0,n1+n2+1])
%     xlabel('Time from others'' landings (sec)')
%     set(gca,'YTick',[],'TickDir','out','YColor','none','FontSize',14)
% 
%     % the arrows
%     xO = -0.2;  
%     yO = 0.5;
%     x_fd1 = t_margin(1)-0.1;
%     x_fd2 = x_fd1;
%     y_fd1 = round(n2+n1/2);
%     y_fd2 = round(n2/2);
%     % patch(...
%     %     [x_fd1+xO x_fd2+xO; x_fd1+xO x_fd1+xO; x_fd1 x_fd2], ...
%     %     [y_fd1-yO y_fd2-yO; y_fd1+yO y_fd2+yO; y_fd1 y_fd2], 'k', 'clipping', 'off')
% 
%     % Texts
%     text(x_fd1-0.1, y_fd1, 'Social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
%     text(x_fd2-0.1, y_fd2, 'Non-social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
% 
%     %%% Save figure
%     saveDir = fullfile(rootdir,'analysis','figure','social_modulation_all', ...
%         sprintf('ephys%s',bname_ephys),trjName,'others_landing',region);
%     if isempty(dir(saveDir)); mkdir(saveDir); end
%     saveas(fig,fullfile(saveDir,sprintf('Unit%d.jpg',unit_id)))
%     close(fig)
% end   
% 
% saveFile = fullfile(rootdir,'analysis','social_modulation.mat');
% if isempty(dir(saveFile))
%     save(saveFile,'mod_other_land');
% else
%     save(saveFile,'mod_other_land','-append');
% end