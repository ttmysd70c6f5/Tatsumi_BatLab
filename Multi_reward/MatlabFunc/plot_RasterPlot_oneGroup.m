function plot_RasterPlot_oneGroup(sptimes_1)

n_trials_1 = length(sptimes_1); % number of trials
% n_trials_2 = length(sptimes_2); % number of trials

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
% for ll = 1:n_trials_2
%     sptimes_trial = sptimes_2(ll).sptimes; % spike times in the current trial
%     % splocs_trial = repmat(ll,size(sptimes_trial)); % location of spikes
%     for ss = 1:length(sptimes_trial)
%         line([sptimes_trial(ss) sptimes_trial(ss)], [cnt-1 cnt],'Color','k')
%         % plot(sptimes_trial,splocs_trial,'.k','MarkerSize', 5);
%     end
%     cnt = cnt + 1;
% end

% yregion(0,n_trials_1+n_trials_2,'FaceColor',[100 100 200]/255)
% yregion(0,n_trials_1,'FaceColor',[200 100 100]/255)

hold off
