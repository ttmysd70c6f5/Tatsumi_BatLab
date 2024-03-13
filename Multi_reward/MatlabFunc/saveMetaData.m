function saveMetaData()
ParentDir = 'Z:\users\Tatsumi\data\MultiReward'; % the directory containing all datasets
recListPath = 'G:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info

allRecs = GetRecList(recListPath); % load the rec info list
allRecNames = unique(allRecs.recName);

for rr = 1:length(allRecNames)
    recname = allRecNames{rr};
    recID = find(strcmp(allRecs.recName,recname),1,'first');

    % extract meta data from the sheet
    clear mtData
    % add basic info
    mtData.recDate = allRecs.date(recID); % recording date
    mtData.rootdir = allRecs.rootdir{recID}; % root directory
    mtData.recName = allRecs.recName{recID}; % name of recording
    % add behavioral settings
    mtData.expType = allRecs.expType{recID}; % type of experiment: social/single bat
    mtData.rewardType = allRecs.rewardType{recID}; % type of reward: single/multi flavor
    mtData.reward_list = strsplit(allRecs.reward_list{recID},','); % reward list
    id_s1 = find(strcmp(allRecs.Properties.VariableNames,'s1'));
    n_session = 4; % initialize the number of session
    for i = 1:4 % s1 ~ s4
        mtData.session(i).reward_session = cellfun(@str2num,strsplit(allRecs{recID,id_s1+i-1}{:},','));
        if sum(mtData.session(i).reward_session) == 0
            n_session = min(n_session,i); % get the first session without feeder activation
        end
    end
    mtData.n_session = n_session - 1; % number of session
    % id_start = find(strcmp(allRecs.Properties.VariableNames,'14447'));
    % id_end = find(strcmp(allRecs.Properties.VariableNames,'9776'));
    bnames = allRecs.Properties.VariableNames(12:21);
    bIncluded = logical(table2array(allRecs(recID,12:21)));
    bnames = bnames(bIncluded);
    mtData.bname = bnames; % bat names
    % add ehys settings
    % search for all the entries within the same recording
    recIDs = find(strcmp(allRecs.rootdir,mtData.rootdir)); % index of the same recording
    for id = 1:length(recIDs)
        recID_now = recIDs(id);
        mtData.ephys(id).ephysDir = allRecs.ephysDir{recID_now}; % ephys directory
        mtData.ephys(id).recType = allRecs.recType{recID_now}; % type of recording: tethered/untethered
        mtData.ephys(id).bname_ephys = allRecs.ephysBat{recID_now}; % name of ephys bat
        mtData.ephys(id).probeNum = allRecs.probeNum(recID_now); % probe number 
        mtData.ephys(id).trjName = allRecs.trjName{recID_now}; % trajectory name
        mtData.ephys(id).hemisphere = allRecs.hemisphere{recID_now}; % trajectory name
    end
    
    
    % save
    if isempty(dir(fullfile(ParentDir,mtData.rootdir,'analysis'))); mkdir(fullfile(ParentDir,mtData.rootdir,'analysis')); end
    save(fullfile(ParentDir,mtData.rootdir,'analysis',sprintf('metadata_%s.mat',mtData.recName)),'mtData')
end