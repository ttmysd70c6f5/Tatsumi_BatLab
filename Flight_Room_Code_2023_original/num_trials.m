% function summary = num_trials(bhv_data)

trial = bhv_data.trials;
field = bhv_data.fields;
bid = trial(:,9);
feeder = trial(:,10);

n_bat = 10;
n_feeder = 4;

summary = zeros(n_bat,n_feeder);
for bb = 1:n_bat
    for ff = 1:n_feeder
        summary(bb,ff) = sum(bid == bb & feeder == ff);
    end
end