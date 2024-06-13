function plot_firing_modulation_twoGroup(t_onset_sec_1,t_end_sec_1,t_onset_sec_2,t_end_sec_2,t_margin,gname,sptimes_unit_sec)
n_flights_1 = length(t_onset_sec_1);
n_flights_2 = length(t_onset_sec_2);

% Sort flights according to the flight length
t_diff_1 = t_end_sec_1 - t_onset_sec_1;
t_diff_2 = t_end_sec_2 - t_onset_sec_2;

[~,sort_idx_1] = sort(t_diff_1,'ascend');
[~,sort_idx_2] = sort(t_diff_2,'ascend');
t_onset_sec_1 = t_onset_sec_1(sort_idx_1);
t_onset_sec_2 = t_onset_sec_2(sort_idx_2);
t_end_sec_1 = t_end_sec_1(sort_idx_1);
t_end_sec_2 = t_end_sec_2(sort_idx_2);

% Extract spike times
[sptimes_sec_1,t_event_sec_1] = alignSpikeTimes_OneUnit(sptimes_unit_sec,t_onset_sec_1,t_end_sec_1,t_margin);
[sptimes_sec_2,t_event_sec_2] = alignSpikeTimes_OneUnit(sptimes_unit_sec,t_onset_sec_2,t_end_sec_2,t_margin);

%%% FIGURE 1: raster plot around landings
nexttile

hold on
plot_RasterPlot_twoGroup(sptimes_sec_1,sptimes_sec_2)

for ff = 1:n_flights_1
    line([0 0], [ff-1 ff],'Color',[250 50 50 255]/255,'LineWidth',2) % onset of trial
    line([t_event_sec_1(ff,2) t_event_sec_1(ff,2)],[ff-1 ff],'Color',[150 150 150 255]/255,'LineWidth',2) % end of trial
end
for ff = 1:n_flights_2
    line([0 0], [n_flights_1+ff-1 n_flights_1+ff],'Color',[150 150 150 255]/255,'LineWidth',2) % onset of trial
    line([t_event_sec_2(ff,2) t_event_sec_2(ff,2)],[n_flights_1+ff-1 n_flights_1+ff],'Color',[150 150 150 255]/255,'LineWidth',2) % end of trial
end
hold off

xlim([t_margin(1), max(max(t_diff_1),max(t_diff_2))+t_margin(2)])
            
xlabel('Time (sec)')
ylabel('Flights')
set(gca,'FontSize',14)

%%% FIGURE 2: Plot firing rate curve
nexttile

sp_bin = 0.1; % bin size (sec)
sp_edges = t_margin(1)-sp_bin/2:sp_bin:t_margin(2)+sp_bin/2; % bin sptimes into 0.1ms bins
% if ~isempty(t_base)
%     t_baseline = t_base(1):sp_bin:t_base(2); % base line to normalize firing rates
% end
plot_firingCurve_TwoGroup(sptimes_sec_1,sptimes_sec_2,sp_edges,sp_bin,gname)

xlabel('Time (sec)')
ylabel('Firing rate')
set(gca,'FontSize',14)
xlim([t_margin(1), max(max(t_diff_1),max(t_diff_2))+t_margin(2)])

