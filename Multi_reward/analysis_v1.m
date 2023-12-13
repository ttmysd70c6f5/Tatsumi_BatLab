% directory
datadir = '/Users/tatsumi/Documents/Yartsev_lab/Data/MultiReward/group_a/recording/231204';

reward_dir  =   fullfile(datadir,'reward');
reward_file =   dir(fullfile(reward_dir,'*.mat'));
cortex_dir  =   fullfile(datadir,'cortex');
cdp_dir     =   fullfile(datadir,'cdp');

% Serial number
s_numbers   =   load(fullfile(datadir,'Serial_numbers_in_use.txt'));
n_bats      =   length(s_numbers);

% Reward timings
n_session =   length(reward_file);
n_feeder  =   4;
raw_summary =   struct([]);
for dd = 1:n_session
    load(fullfile(reward_file(dd).folder,reward_file(dd).name),'bhv_data');
    raw_summary_now =   zeros(n_bats,n_feeder);
    bid     =   bhv_data.trials(:,9);
    bfeeder =   bhv_data.trials(:,10);
    for bb = 1:n_bats
        for ff = 1:n_feeder
            raw_summary_now(bb,ff)  =   sum(bid==bb&bfeeder==ff);
        end
    end
    raw_summary(dd).reward  =   raw_summary_now;
end