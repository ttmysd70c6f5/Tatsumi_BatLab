function varargout = LoadExperimentData(recID,dataTypes)

% Specify the index of session to be loaded
%% 240207
% recID = 3;
% recID = 4;
% recID = 5;
%% 240208
% recID = 6;
% recID = 7;
% recID = 8;
%% 240209_1
% recID = 9;
% recID = 10;
%% 240209_2
% recID = 11;
%% 240210_1
% recID = 12;
% recID = 13;
%% 240210_2
% recID = 14;
%% 240211_1
% recID = 15;
% recID = 16;
% recID = 17;
%% 240211_2
% recID = 18;
%% 240212
% recID = 19;
% recID = 20;
% recID = 21;

% Load the table of rec info
ParDir = 'Z:\users\Tatsumi\data\MultiReward'; % parent directory
recListPath = 'H:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info
allRecs = GetRecList(recListPath); % load the rec info list

% rootdir
rootdir = fullfile(ParDir,allRecs.rootdir{recID});
fprintf('Loading data from %s ...\n',rootdir)
cd(rootdir)

% ephys directory
if ~isempty(allRecs.ephysDir)
    ephysDir = fullfile(rootdir,'ephys',allRecs.ephysDir{recID},allRecs.ephysBat{recID},'Processed');
    probeNum = allRecs.probeNum(recID); % probe number
end

% Meta data
expType = allRecs.expType{recID}; % experiment type: social/solo/tethered
rewardType = allRecs.rewardType{recID}; % reward type: multi/single flavor
rootdir = fullfile(ParDir,allRecs.rootdir{recID}); % root directory of the recording
% ephysDir = fullfile(rootdir,'ephys',allRecs.ephysDir{recID}); % ephys directory
recDate = char(allRecs.date(recID),'MM/dd/yy'); % date of recording
recName = allRecs.recName{recID};
bname_ephys = allRecs.ephysBat{recID}; % name of ephys bat
probeNum = allRecs.probeNum(recID); % probe number
trjName = allRecs.trjName{recID}; % trajectory

% Create output directories
outDir = fullfile(rootdir,'analysis'); % output
figDir = fullfile(rootdir,'analysis','figure'); % figure
if isempty(dir(outDir)); mkdir(outDir); end
if isempty(dir(figDir)); mkdir(figDir); end

% Adhoc brain region identification
Layers = struct([]);

%%% 14633, MEC_2
Layers(1).batid = '14633';
Layers(1).trjName = 'MEC_2'; % name of trajectory

Layers(1).layer(1).name = 'MEC';
Layers(1).layer(1).depth = [0,2800];
Layers(1).layer(1).color = [127,201,127]/255;

Layers(1).layer(2).name = 'VC';
Layers(1).layer(2).depth =[5500,8000];
Layers(1).layer(2).color = [190,174,212]/255;

Layers(1).layer(3).name = 'unknown';
Layers(1).layer(3).depth = [0,0];
Layers(1).layer(3).color = [0,0,0];

%%% 14662, dCA1_MEC_3
Layers(2).batid = '14662'; % id of bat
Layers(2).trjName = 'dCA1_MEC_3'; % name of trajectory
Layers(2).layer(1).name = 'MEC';
Layers(2).layer(1).depth = [0,500];
Layers(2).layer(1).color = [127,201,127]/255;

Layers(2).layer(2).name = 'Subiculum';
Layers(2).layer(2).depth =[500,1600];
Layers(2).layer(2).color = [190,174,212]/255;

Layers(2).layer(3).name = 'DG';
Layers(2).layer(3).depth =[1600,4100];
Layers(2).layer(3).color = [253,192,134]/255;

Layers(2).layer(4).name = 'Subiculum';
Layers(2).layer(4).depth =[4100,5200];
Layers(2).layer(4).color = [190,174,212]/255;

Layers(2).layer(5).name = 'VC';
Layers(2).layer(5).depth =[4100,5200];
Layers(2).layer(5).color = [150,150,150]/255;

Layers(2).layer(6).name = 'unknown';
Layers(2).layer(6).depth = [0,0];
Layers(2).layer(6).color = [0,0,0];


%%% 14662, dCA1_vCA1_3
Layers(3).batid = '14662'; % id of bat
Layers(3).trjName = 'dCA1_vCA1_3'; % name of trajectory

Layers(3).layer(1).name = 'CA1';
Layers(3).layer(1).depth =[0,1700];
Layers(3).layer(1).color = [190,174,212]/255;

Layers(3).layer(2).name = 'CA3';
Layers(3).layer(2).depth = [1700, 3300];
Layers(3).layer(2).color = [253,192,134]/255;

Layers(3).layer(3).name = 'CA1';
Layers(3).layer(3).depth =[3300,4100];
Layers(3).layer(3).color = [190,174,212]/255;

Layers(3).layer(4).name = 'unknown';
Layers(3).layer(4).depth = [0,0];
Layers(3).layer(4).color = [0,0,0];



% load data
n = length(dataTypes);
% i = 0; % current argout size
for k = 1:n % seperate the simbol for data types
    str = dataTypes(k); % current data type
    switch str
        case 'a' % load analog signals
            disp('Loading behavioral data ...')
            [~, varargout{k}] = PreprocessBhvData(rootdir);
        case 'b' % load behaviors
            disp('Loading analog signals ...')
            [varargout{k}, ~] = PreprocessBhvData(rootdir);
        case 'c'
        case 'l'
            disp('Loading LFP data ...')
            varargout{k} = loadTrodesLFP(ephysDir);
        case 'm'
            disp('Loading meta data ...')
            metaDir = dir(fullfile(rootdir,'analysis','metadata*.mat'));
            load(fullfile(metaDir.folder,metaDir.name),'mtData');
            varargout{k} = mtData;
            mtData.bname_ephys = allRecs.ephysBat{recID}; % name of ephys bat
            mtData.trjName = allRecs.trjName{recID}; % trajectory
        case 's'
            disp('Loading spike data ...')
            spikeData = loadSpikeData(ephysDir,'kilosort_out_folder_name', sprintf('kilosort_outdir_probe%d',probeNum));            
            
            % Assign region to the spike data structure
            unit_ids = spikeData.good_units.cluster_id;
            regions = cell(length(unit_ids),1);
            for unit = 1:length(unit_ids)
                depth = spikeData.good_units.depth(unit); % depth of unit
                layers = Layers(strcmp({Layers.trjName},trjName) & strcmp({Layers.batid},bname_ephys)); % specify the layer
                region = get_probe_region(layers,depth);
                regions{unit} = region;
            end
            spikeData.good_units.region = regions;
            varargout{k} = spikeData;
        case 'r' % load reward data
            disp('Loading reward data ...')
            varargout{k} = PreprocessRewardData(rootdir);
    end
end

% Landings = extractLandings(BhvData,RewardData,mtData,'include','all');