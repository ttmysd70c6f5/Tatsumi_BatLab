function plot_firingCurve_oneGroup(sptimes_sec_1,sp_edges,sp_bin)
% Parameters
filt_smp = 5; % number of samples to smooth

% Binning the spike times
% [fr_1,sem_1] = compute_psth(sptimes_sec_1,sp_edges,sp_bin,filt_smp);
[fr_1,sem_1] = compute_psth_v2(sptimes_sec_1,sp_edges,sp_bin,filt_smp);
% [fr_2,sem_2] = compute_psth(sptimes_sec_2,sp_edges,filt_smp);
% N_1 = zeros(length(sptimes_sec_1),length(sp_edges)-1);
% N_2 = zeros(length(sptimes_sec_2),length(sp_edges)-1);
% 
% % group 1
% for ff = 1:length(sptimes_sec_1)
%     [N,~] = histcounts(sptimes_sec_1(ff).sptimes,sp_edges);
%     N_1(ff,:) = N;
% end
% N_1 = smoothdata(N_1,2,'gaussian',filt_smp); % gaussian filtering
% sem_1 = std(N_1,0,1) / sqrt(size(N_1,1)); % standard error of mean
% fr_1 = mean(N_1,1); % 
% 
% % group 2
% for ff = 1:length(sptimes_sec_2)
%     [N,~] = histcounts(sptimes_sec_2(ff).sptimes,sp_edges);
%     N_2(ff,:) = N;
% end
% N_2 = smoothdata(N_2,2,'gaussian',filt_smp);
% sem_2 = std(N_2,0,1,'includemissing') / sqrt(size(N_2,1));
% fr_2 = mean(N_2,1);

% y_max = max(max(fr_1+sem_1),max(fr_2 + sem_2)); 
% y_min = min(min(fr_1-sem_1),min(fr_2 - sem_2));
y_max = max(fr_1+sem_1); 
y_min = min(fr_1-sem_1);

%==== figure
hold on
% group 1
x = mean([sp_edges(1:end-1); sp_edges(2:end)],1);
% x = sp_edges(1)+sp_bin/2:sp_bin:sp_edges(end)-sp_bin/2;
p1 = boundedline(x,fr_1,sem_1,'r','alpha','transparency', 0.1);
% group 2
% p2 = boundedline(x,fr_2,sem_2,'b','alpha','transparency', 0.1);


% line at t = 0
% p2 = line([0 0], [0 y_max],'Color',[150 150 150 255]/255,'LineWidth',2);
hold off

% legend(p1,gname)

if y_max ~= y_min
    % ylim([(3*y_min-y_max)/2 (3*y_max-y_min)/2])
    ylim([0,y_max*1.3])
end
