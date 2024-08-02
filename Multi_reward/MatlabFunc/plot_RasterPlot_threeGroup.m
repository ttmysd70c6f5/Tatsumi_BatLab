function plot_RasterPlot_threeGroup(sptimes_1,sptimes_2,sptimes_3,g_clr)

n_trials_1 = length(sptimes_1); % number of trials
n_trials_2 = length(sptimes_2); % number of trials
n_trials_3 = length(sptimes_3); % number of trials

% sort trials according to the event length
% [~,idx] = sort(sort_vec,'ascend'); % sort trials
% t_onset_1 = t_onset_1(idx);
% t_end_1 = t_end_1(idx);
% sptimes_1 = sptimes_1(idx);

% figure
cnt = 1; % count
hold on
for ll = 1:n_trials_1
    sptimes_trial = sptimes_1(ll).sptimes; % spike times in the current trial
    % splocs_trial = repmat(ll,size(sptimes_trial)); % location of spikes
    for ss = 1:length(sptimes_trial)
        line([sptimes_trial(ss) sptimes_trial(ss)], [cnt-1 cnt],'Color','k')
        % plot(sptimes_trial,splocs_trial,'.k','MarkerSize', 5);
    end
    cnt = cnt + 1;
end
for ll = 1:n_trials_2
    sptimes_trial = sptimes_2(ll).sptimes; % spike times in the current trial
    % splocs_trial = repmat(ll,size(sptimes_trial)); % location of spikes
    for ss = 1:length(sptimes_trial)
        line([sptimes_trial(ss) sptimes_trial(ss)], [cnt-1 cnt],'Color','k')
        % plot(sptimes_trial,splocs_trial,'.k','MarkerSize', 5);
    end
    cnt = cnt + 1;
end
for ll = 1:n_trials_3
    sptimes_trial = sptimes_3(ll).sptimes; % spike times in the current trial
    % splocs_trial = repmat(ll,size(sptimes_trial)); % location of spikes
    for ss = 1:length(sptimes_trial)
        line([sptimes_trial(ss) sptimes_trial(ss)], [cnt-1 cnt],'Color','k')
        % plot(sptimes_trial,splocs_trial,'.k','MarkerSize', 5);
    end
    cnt = cnt + 1;
end

yregion(0,n_trials_1+n_trials_2+n_trials_3,'FaceColor',g_clr(3,:)/255)
yregion(0,n_trials_1+n_trials_2,'FaceColor',g_clr(2,:)/255)
yregion(0,n_trials_1,'FaceColor',g_clr(1,:)/255)

hold off