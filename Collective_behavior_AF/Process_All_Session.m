function Process_All_Session()

%% DBSCAN
%% Clustering resting positions for all sessions
% options
opt.smpts = 1000000; % sampling size
opt.minpts = 100; % parameter of dbscan
opt.th = 0.01; % threshold for exclusion of small clusters

% DBSCAN with knn-search
tic
DBSCAN_knnsearch_clustering(opt,false)
toc


%% State transitions analysis
state_transitions_sessions()

end