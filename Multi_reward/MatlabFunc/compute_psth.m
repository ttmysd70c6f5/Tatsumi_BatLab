function [fr,sem] = compute_psth(sptimes,sp_edges,sp_bin,filt_smp)
% Binning the spike times
N = zeros(length(sptimes),length(sp_edges)-1);

% group 1
for ff = 1:length(sptimes)
    [N(ff,:),~] = histcounts(sptimes(ff).sptimes,sp_edges);
end
N = smoothdata(N,2,'movmean',filt_smp) / sp_bin;
fr = mean(N,1);
sem = std(N,0,1) / sqrt(size(N,1));

% N = N/sp_bin;
% fr = mean(N,1);
% sem = std(N,0,1) / sqrt(size(N,1)); % standard error of mean
% fr = smoothdata(fr,2,'gaussian',filt_smp); % gaussian filtering

% if ~isempty(t_baseline)
%     bin_normalize = ismember(round(mean([sp_edges(1:end-1);sp_edges(2:end)],1),2), round(t_baseline,2));
%     mean_base = mean(fr(bin_normalize));
%     std_base = std(fr(bin_normalize),0,2);
%     fr = (fr - mean_base) / std_base;
% end