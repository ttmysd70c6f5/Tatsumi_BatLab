function varargout = loadExpData(ParentDir,allRecs,recID,dataTypes)
% loadExpData(AllRecTable,dataTypes) 
% Input
% ParentDir: the parent directory for all recordings
% AllRecTable: the table of the recording infor
% recID: Row number of the data to be loaded
% dataType: the type of data to load, specified with the concatenated simbols.
%           a = Analog signals, b = behaviors, c = camera, l = lfp, m = metadata, 
%           r = reward, s = spikes

% Input variables
p = inputParser;
addRequired(p,'ParentDir');
addRequired(p,'allRecs');
addRequired(p,'recID');
addRequired(p,'dataType');

parse(p,ParentDir,allRecs,recID,dataTypes)


% rootdir
rootdir = fullfile(ParentDir,allRecs.rootdir{recID});
fprintf('Loading data from %s ...\n',rootdir)

% ephys directory
if ~isempty(allRecs.ephysDir)
    ephysDir = fullfile(rootdir,'ephys',allRecs.ephysDir{recID},allRecs.ephysBat{recID},'Processed');
    probeNum = allRecs.probeNum(recID); % probe number
end

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
        case 's'
            disp('Loading spike data ...')
            varargout{k} = loadSpikeData(ephysDir,'kilosort_out_folder_name', sprintf('kilosort_outdir_probe%d',probeNum));

        case 'r' % load reward data
            disp('Loading reward data ...')
            varargout{k} = PreprocessRewardData(rootdir);
    end
end

end