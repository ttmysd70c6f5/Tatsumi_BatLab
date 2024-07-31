%% Load neural data
recListPath = 'H:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info
allRecs = GetRecList(recListPath); % load the rec info list

% Load data
% recID_1 = 3; recID_2 = 4; % 240207, #14633 and 14462, MEC
% recID_1 = 6; recID_2 = 7; % 240208, #14633 and 14462, MEC
% recID_1 = 3; recID_2 = 4; % 240209, #14633 and 14462, MEC
% recID_1 = 12; recID_2 = 13; % 240210_1, #14633 and 14462, MEC
% recID_1 = 3; recID_2 = 4; % 240212, #14633 and 14462, MEC
% dataTypes = 'barslm'; % data types to be loaded
% [BhvData_1,AnalogData_1,RewardData_1,spikeData_1,lfpData_1,mtData_1] = LoadExperimentData(recID_1,dataTypes);
% [BhvData_2,AnalogData_2,RewardData_2,spikeData_2,lfpData_2,mtData_2] = LoadExperimentData(recID_2,dataTypes);
recID = 12;
dataTypes = 'barsm'; % data types to be loaded
[BhvData,AnalogData,RewardData,spikeData,mtData] = LoadExperimentData(recID,dataTypes);
[Flights, TrjCluster] = extractFlights(BhvData,RewardData,mtData,'include','all');
%%
bname_ephys = allRecs.ephysBat{recID}; % name of ephys bat
mtData.bname_ephys = bname_ephys;
mtData.trjName = allRecs.trjName{recID}; % trajectory
observe_bat = find(strcmp(bname_ephys,mtData.bname)); % Ephys bat
% observe_bat = 6;

%% Plot firing around events
%%% Plot firing around self takeoffs
% Directory to save figures
target = 'takeoff';
FigDir = fullfile('F:\Data\GroupForaging\Analysis\figure',sprintf('Social_modulation_self_%s',target));

% Params
t_margin = [-1 1];
fly_bat = find(strcmp(mtData.bname_ephys,mtData.bname)); % Ephys bat

% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0 0.05]);
PlotClusteredTrajectory(TrjCluster,mtData,fullfile(FigDir,'trj_cluster'))

% Exclude flights outside ephys recording
t_ephys_min = min(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(1);
t_ephys_max = max(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(2);
t_plot_min = BhvData.t(1)-t_margin(1);
t_plot_max = BhvData.t(end)-t_margin(2);

SocialFlights = ExtractSelfSocialFlights_v3(BhvData,Flights,TrjCluster,fly_bat,target,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));

saveDir = fullfile(FigDir,mtData.recName,sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
PlotSelfFlightModulation(BhvData,SocialFlights,spikeData,target,saveDir)
PlotFlightModulation_notCluster(BhvData,SocialFlights,spikeData,target,saveDir)

load(fullfile(saveDir,'pefr_all.mat'))
fr_all = [fr_all, repmat({'self_toff',mtData.recName},[size(fr_all,1),1])];
fr_field = [fr_field, {'type','recname'}];
self_toff_pefr_all = fr_all;
save(fullfile(saveDir,'pefr_all.mat'),'fr_all','fr_field')

%%

%%% Plot firing around self landings
% Directory to save figures
target = 'landing';
FigDir = fullfile('F:\Data\GroupForaging\Analysis\figure',sprintf('Social_modulation_self_%s',target));

% Params
t_margin = [-1 1];
fly_bat = find(strcmp(mtData.bname_ephys,mtData.bname)); % Ephys bat

% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0 0.05]);
PlotClusteredTrajectory(TrjCluster,mtData,fullfile(FigDir,'trj_cluster'))

% Exclude flights outside ephys recording
t_ephys_min = min(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(1);
t_ephys_max = max(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(2);
t_plot_min = BhvData.t(1)-t_margin(1);
t_plot_max = BhvData.t(end)-t_margin(2);

% Extract social flights
SocialFlights = ExtractSelfSocialFlights_v3(BhvData,Flights,TrjCluster,fly_bat,target,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));

% Plot firing and raster plots
saveDir = fullfile(FigDir,mtData.recName,sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
PlotSelfFlightModulation(BhvData,SocialFlights,spikeData,target,saveDir)
PlotFlightModulation_notCluster(BhvData,SocialFlights,spikeData,target,saveDir)

% Load extracted data for t-SNE
load(fullfile(saveDir,'pefr_all.mat'))
fr_all = [fr_all, repmat({'self_land',mtData.recName},[size(fr_all,1),1])];
fr_field = [fr_field, {'type','recname'}];
self_land_pefr_all = fr_all;
save(fullfile(saveDir,'pefr_all.mat'),'fr_all','fr_field')

%%
%%% Plot firing around others' takeoffs
% Directory to save figures
target = 'takeoff';
FigDir = fullfile('F:\Data\GroupForaging\Analysis\figure',sprintf('Social_modulation_other_%s',target));

% Params
t_margin = [-1 1];

% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0 0.05]);
PlotClusteredTrajectory(TrjCluster,mtData,fullfile(FigDir,'trj_cluster'))

% [SocialFlights, AsocialFlights] = ExtractSocialFlights(BhvData,Flights,observe_bat,target);
% [SocialFlights, AsocialFlights] = ExtractSocialFlights_v2(BhvData,Flights,TrjCluster,observe_bat,target);

% Exclude flights outside ephys recording
t_ephys_min = min(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(1);
t_ephys_max = max(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(2);
t_plot_min = BhvData.t(1)-t_margin(1);
t_plot_max = BhvData.t(end)-t_margin(2);

SocialFlights = ExtractOthersSocialFlights_v3(BhvData,Flights,TrjCluster,observe_bat,target,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));
% SocialFlights = CorrectFlights(SocialFlights,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));

saveDir = fullfile(FigDir,mtData.recName,sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
PlotOtherFlightModulation(BhvData,SocialFlights,spikeData,target,saveDir);
PlotFlightModulation_notCluster(BhvData,SocialFlights,spikeData,target,saveDir)

load(fullfile(saveDir,'pefr_all.mat'))
fr_all = [fr_all, repmat({'oth_toff',mtData.recName},[size(fr_all,1),1])];
fr_field = [fr_field, {'type','recname'}];
oth_toff_pefr_all = fr_all;
save(fullfile(saveDir,'pefr_all.mat'),'fr_all','fr_field')


%%% Plot firing around others' landings
% Directory to save figures
target = 'landing';
FigDir = fullfile('F:\Data\GroupForaging\Analysis\figure',sprintf('Social_modulation_other_%s',target));

% Params
t_margin = [-1 1];

% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0.95 1]);
PlotClusteredTrajectory(TrjCluster,mtData,fullfile(FigDir,'trj_cluster'))

% [SocialFlights, AsocialFlights] = ExtractSocialFlights(BhvData,Flights,observe_bat,target);
% [SocialFlights, AsocialFlights] = ExtractSocialFlights_v2(BhvData,Flights,TrjCluster,observe_bat,target);

% Exclude flights outside ephys recording
t_ephys_min = min(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(1);
t_ephys_max = max(spikeData.allGlobalSpikeTimes_usec/1e6)-t_margin(2);
t_plot_min = BhvData.t(1)-t_margin(1);
t_plot_max = BhvData.t(end)-t_margin(2);

SocialFlights = ExtractOthersSocialFlights_v3(BhvData,Flights,TrjCluster,observe_bat,target,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));
% SocialFlights = CorrectFlights(SocialFlights,max([t_ephys_min,t_plot_min]),min([t_ephys_max,t_plot_max]));

saveDir = fullfile(FigDir,mtData.recName,sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
PlotOtherFlightModulation(BhvData,SocialFlights,spikeData,target,saveDir);
PlotFlightModulation_notCluster(BhvData,SocialFlights,spikeData,target,saveDir)

load(fullfile(saveDir,'pefr_all.mat'))
fr_all = [fr_all, repmat({'oth_land',mtData.recName},[size(fr_all,1),1])];
fr_field = [fr_field, {'type','recname'}];
oth_land_pefr_all = fr_all;
save(fullfile(saveDir,'pefr_all.mat'),'fr_all','fr_field')



%% t-SNE
% Load data
saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\Social_modulation_self_takeoff',sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
load(fullfile(saveDir,'pefr_all.mat'))
self_toff_pefr_all = fr_all;
self_toff_pefr_all = [self_toff_pefr_all, repmat({'self_toff'},[size(self_toff_pefr_all,1),1])];

saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\Social_modulation_self_landing',sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
load(fullfile(saveDir,'pefr_all.mat'))
self_land_pefr_all = fr_all;
self_land_pefr_all = [self_land_pefr_all, repmat({'self_land'},[size(self_land_pefr_all,1),1])];

saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_takeoff',sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
load(fullfile(saveDir,'pefr_all.mat'))
oth_toff_pefr_all = fr_all;
oth_toff_pefr_all = [oth_toff_pefr_all, repmat({'oth_toff'},[size(oth_toff_pefr_all,1),1])];

saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_landing',sprintf('ephys%s_%s',mtData.bname_ephys,mtData.trjName));
load(fullfile(saveDir,'pefr_all.mat'))
oth_land_pefr_all = fr_all;
oth_land_pefr_all = [oth_land_pefr_all, repmat({'oth_land'},[size(oth_land_pefr_all,1),1])];

% Concatenate data
fr_all = [self_toff_pefr_all; self_land_pefr_all; oth_toff_pefr_all; oth_land_pefr_all];

% Exclusion criteria
n_samples = size(fr_all,1);
fr_all = fr_all(cellfun(@(x) ~isempty(x),fr_all(:,7)),:);
idx = cell2mat(fr_all(:,3)) > 1 & cell2mat(fr_all(:,5)) < 0.01/n_samples & cell2mat(fr_all(:,6)) < 0.01/n_samples & cell2mat(fr_all(:,7)) < 0.01/n_samples; % fr, ctest, anova, glm with Bonferroni correction
fprintf('%d/%d cells are responsive.\n',sum(idx),n_samples)
fr_analyze = fr_all(idx,:);

%% Distribution of firing rate
figure
histogram(cell2mat(fr_analyze(:,3)),0:1:50)
xlabel('Firing rate (Hz)')
ylabel('Count')
set(gca,'FontSize',16)

%%
% groups
g_animal = fr_analyze(:,1); % animal id
g_clu = cellfun(@(x) num2str(x),fr_analyze(:,4),'UniformOutput',false); % cluster id
g_unit = cellfun(@(x) num2str(x),fr_analyze(:,2),'UniformOutput',false); % unit id
g_ani_clu = cellfun(@(x,y) join([x y],'_'),g_animal,g_clu,'UniformOutput',false); % animal id and cluster id
g_fr = cell(size(fr_analyze(:,3))); g_fr(:) = {'high'}; g_fr(cell2mat(fr_analyze(:,3))<35) = {'middle-high'}; g_fr(cell2mat(fr_analyze(:,3))<25) = {'middle-low'}; g_fr(cell2mat(fr_analyze(:,3))<15) = {'low'}; 
g_session = fr_analyze(:,9);
g_self = fr_analyze(:,9); 
g_self(strcmp(g_self,'self_toff')) = {'self'}; g_self(strcmp(g_self,'self_land')) = {'self'};
g_self(strcmp(g_self,'oth_toff')) = {'other'}; g_self(strcmp(g_self,'oth_land')) = {'other'};

% Apply t-SNE
X = cell2mat(cellfun(@(x) x',fr_analyze(:,8),'UniformOutput',false));
% Y = tsne(X);
Y = tsne(X,'Algorithm','exact','Distance','euclidean');

% Visualize the result
gr = g_animal;
figure
gscatter(Y(:,1),Y(:,2),gr)
xlabel('Dim 1')
ylabel('Dim 2')

gr = g_unit;
figure
gscatter(Y(:,1),Y(:,2),gr)
xlabel('Dim 1')
ylabel('Dim 2')

gr = g_fr;
figure
gscatter(Y(:,1),Y(:,2),gr)
xlabel('Dim 1')
ylabel('Dim 2')

gr = g_session;
figure
gscatter(Y(:,1),Y(:,2),gr)
xlabel('Dim 1')
ylabel('Dim 2')

gr = g_self;
figure
gscatter(Y(:,1),Y(:,2),gr)
xlabel('Dim 1')
ylabel('Dim 2')

%% Plot firing curve and raster plots for each cluster
% Parameters




% Load behavioral features
% r = BhvData.r(:,:,bb); % position
% v = BhvData.v_abs(:,bb); % velocity
a_abs = BhvData.ac_abs(:,bb); % acceleration
t = BhvData.t; % time
% n_bats = BhvData.n_tags; % # of bats
% blanding = 1 - BhvData.bflying; % 1 if at resting sites
w_step_2 = BhvData.fs; % sampling frequency for behavioral data
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
            ps_ctest = NaN(1,n_unit);
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
                ps_ctest(unit) = stat_test_poisson(3,opt_bin+pts_margin); % p-value of the optimal bin
    
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
                spmat2_social = NaN(n1,diff(t_margin)*w_step_2+1);
                spmat2_asocial = NaN(n2,diff(t_margin)*w_step_2+1);
                for ff = 1:n1
                    t_edges_2 = t_social(ff)+t_margin(1)-1/w_step_2/2:1/w_step_2:t_social(ff)+t_margin(2)+1/w_step_2/2;
                    spmat2_social(ff,:) = histcounts(sptimes_unit_sec,t_edges_2);
                end
                for ff = 1:n2
                    t_edges_2 = t_asocial(ff)+t_margin(1)-1/w_step_2/2:1/w_step_2:t_asocial(ff)+t_margin(2)+1/w_step_2/2;
                    spmat2_asocial(ff,:) = histcounts(sptimes_unit_sec,t_edges_2);
                end
                
                % Peri-event firing rate (social vs non-social)
                pefr = sum(spmat2_social,'all')/(n1*diff(t_margin));
                pefr_asocial = sum(spmat2_asocial,'all')/(n2*diff(t_margin));
    
    
                %%% FIGURE: Firing rate curve and raster plot
                %%% Initialize figure bundle
                fig = figure('visible','off');
                % fig = figure;
                set(fig,'units','normalize','Position',[0.2 0.2 0.3 0.7])
                tl = tiledlayout(2,1,'TileSpacing','tight');
                title_str = sprintf('Unit%d, %s, %dum, %.1fHz (Social: %.1fHz, Non-social: %.1fHz)',unit_id,region,depth,fr,pefr,pefr_asocial);
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
                t_bins = t_margin(1):1/w_step_2:t_margin(2);
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
            mod_other_toff(bb).p_ctest(clu).p_ctest = ps_ctest;
            mod_other_toff(bb).p_glm(clu).p_glm = p_glm;
        end
    end
    mod_other_toff(bb).gen_fig = gen_fig;
% end
saveFile = fullfile(FigDir,'social_modulation_other_toff.mat');
save(saveFile,'mod_other_toff');


%% Plot motor related activity (velocity, acceleration) distribution
observe_bat = find(strcmp(mtData.bname_ephys,mtData.bname)); % observor

% acceleration distribution
saveDir = fullfile(FigDir,'ac_adjust');
savePrefix = 'ac_before_adjust_';
plotACdistribution(SocialFlights,observe_bat,saveDir,savePrefix)
    
% velocity distribution
saveDir = fullfile(FigDir,'v_adjust');
savePrefix = 'v_before_adjust_';
plotVelocityDistribution(SocialFlights,observe_bat,saveDir,savePrefix)

%% 

%% Adjust the acceleration distribution
bin = 0.25;
seed = 1234;

SocialFlights = AdjustACdistribution(SocialFlights,observe_bat,bin,seed);
%% Plot after adjustment
saveDir = fullfile(FigDir,'ac_adjust');
savePrefix = 'ac_after_adjust_';

% acceleration distribution
saveDir = fullfile(FigDir,'ac_adjust');
savePrefix = 'ac_after_adjust_';
plotAdjustedACdistribution(SocialFlights,observe_bat,saveDir,savePrefix)

%%
% Adjusting accleration distributions and generate firing curve plots and raster plots
genPlot_social_mod_oth_takeoff(BhvData,Flights,TrjCluster,spikeData,mtData,FigDir)