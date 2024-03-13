function AllRecTable = GetRecList(recListPath)
% AllRecTable = GetRecList() returns the list of the all recordings.


% recListPath = 'G:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info
opt=detectImportOptions(recListPath,'ReadVariableNames',true); % get the default import options
opt = setvartype(opt,'recName','char'); % set the variable session as char
opt = setvartype(opt,'ephysBat','char'); % set the variable session as char
opt.PreserveVariableNames = 1; % preserve the variable names
AllRecTable = readtable(recListPath,opt); % load rec info
end