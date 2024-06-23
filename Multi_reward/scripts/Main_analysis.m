%% Load data
clear all
% Load the table of rec info
ParDir = 'Z:\users\Tatsumi\data\MultiReward'; % parent directory
recListPath = 'K:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info
allRecs = GetRecList(recListPath); % load the rec info list

% Specify the index of session to be loaded
%%% 240207
% recID = 3;
% recID = 4;
% recID = 5;
%%% 240208
% recID = 6;
% recID = 7;
% recID = 8;
%%% 240209_1
% recID = 9;
% recID = 10;
%%% 240209_2
% recID = 11;
% %% 240210_1
recID = 12;
% recID = 13;
%%% 240210_2
% recID = 14;
%%% 240211_1
% recID = 15;
% recID = 16;
% recID = 17;
%%% 240211_2
% recID = 18;
%%% 240212
% recID = 19;
% recID = 20;
% recID = 21;


% general information
expType = allRecs.expType{recID}; % experiment type: social/solo/tethered
rewardType = allRecs.rewardType{recID}; % reward type: multi/single flavor
rootdir = fullfile(ParDir,allRecs.rootdir{recID}); % root directory of the recording
ephysDir = fullfile(rootdir,'ephys',allRecs.ephysDir{recID}); % ephys directory
recDate = char(allRecs.date(recID),'MM/dd/yy'); % date of recording
recName = allRecs.recName{recID};
bname_ephys = allRecs.ephysBat{recID}; % name of ephys bat
probeNum = allRecs.probeNum(recID); % probe number
trjName = allRecs.trjName{recID}; % trajectory

cd(rootdir)

% Create output directories
outDir = fullfile(rootdir,'analysis'); % output
figDir = fullfile(rootdir,'analysis','figure'); % figure
if isempty(dir(outDir)); mkdir(outDir); end
if isempty(dir(figDir)); mkdir(figDir); end

% load data
% dataTypes = 'barslm'; % data types to be loaded
% [BhvData,AnalogData,RewardData,spikeData,lfpData,mtData] = loadExpData(ParDir, allRecs, recID, dataTypes); % load behavior, analog signals, and reward data
dataTypes = 'barsm'; % data types to be loaded
[BhvData,AnalogData,RewardData,spikeData,mtData] = loadExpData(ParDir, allRecs, recID, dataTypes); % load behavior, analog signals, and reward data
% dataTypes = 'barm'; % data types to be loaded
% [BhvData,AnalogData,RewardData,mtData] = loadExpData(ParDir, allRecs, recID, dataTypes); % load behavior, analog signals, and reward data
% dataTypes = 's'; % data types to be loaded
% [spikeData] = loadExpData(ParDir, allRecs, recID, dataTypes); % load behavior, analog signals, and reward data
% Landings = extractLandings(BhvData,RewardData,mtData,'include','all');
[Flights, ~] = extractFlights(BhvData,RewardData,mtData,'include','all','start',2*60);

% Meta Data
mtData.bname_ephys = bname_ephys;
mtData.recName = recName;
mtData.trjName = trjName;

%Layer information
spikeData = adhoc_ch_reconstruct(spikeData,trjName,bname_ephys);

%% Agent-vector coding
saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\vector_coding',mtData.recName,sprintf('%s_%s',mtData.bname_ephys,mtData.trjName));
if isempty(dir(saveDir)); mkdir(saveDir); end
social_vector_coding(BhvData,Flights,spikeData,mtData,saveDir)
plot_social_vector_coding(BhvData,spikeData,mtData,saveDir)

%% Agent-vector coding for each landing
% Clustering of trajectories
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0.95 1]);
saveDir = fullfile('F:\Data\GroupForaging\Analysis\figure\vector_coding',mtData.recName,sprintf('%s_%s',mtData.bname_ephys,mtData.trjName),'every_cluster');
FigDir = fullfile(saveDir,'trjcluster');
PlotClusteredTrajectory(TrjCluster,mtData,FigDir)
% Compute the firing map for each cluster
social_vector_coding_every_land(BhvData,Flights,TrjCluster,spikeData,mtData,saveDir)

% Plot firing map
plot_social_vector_coding_every_land(BhvData,TrjCluster,spikeData,mtData,saveDir)

%% Social modulation by others' landings
% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0.95 1]);
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_land\trj_cluster';
PlotClusteredTrajectory(TrjCluster,mtData,FigDir)

% Adjusting accleration distributions and generate firing curve plots and raster plots
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_land\';
genPlot_social_mod_oth_land(BhvData,Flights,TrjCluster,spikeData,mtData,FigDir)

%% Social modulation by others' takeoffs
% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0 0.05]);
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_takeoff\trj_cluster';
PlotClusteredTrajectory(TrjCluster,mtData,FigDir)

% Adjusting accleration distributions and generate firing curve plots and raster plots
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_other_takeoff\';
genPlot_social_mod_oth_takeoff(BhvData,Flights,TrjCluster,spikeData,mtData,FigDir)

%% Social modulation by self landings
% Hierarchical clustering of trajectories based on ending portions
TrjCluster = TrjClustering_HC_Frechet_v2(Flights,'cutoff',1.0,'iprange',[0.95 1]);
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_self_land\trj_cluster';
PlotClusteredTrajectory(TrjCluster,mtData,FigDir)

% Adjusting accleration distributions and generate firing curve plots and raster plots
FigDir = 'F:\Data\GroupForaging\Analysis\figure\Social_modulation_self_land\';
genPlot_social_mod_self_land(BhvData,Flights,TrjCluster,spikeData,mtData,FigDir)