function deviation = compute_fr_deviation_oneUnit(sptimes_sec,fr,t_edge)

% Calculate the # of spikes with a sliding window
spcounts = histcounts(sptimes_sec,t_edge);
deviation = arrayfun(@(x) -log10(poisspdf(x,fr)),spcounts);