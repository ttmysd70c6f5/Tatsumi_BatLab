function deviation = compute_fr_deviation(sptimes_all_usec,fr,t_edge)

n_unit = length(sptimes_all_usec);
deviation = NaN(n_unit,length(t_edge)-1);
for unit = 1:n_unit
    sptimes_sec = sptimes_all_usec{unit}/1e6;
    
    % Calculate the # of spikes with a sliding window
    
    spcounts = histcounts(sptimes_sec,t_edge);
    deviation(unit,:) = arrayfun(@(x) -log10(poisspdf(x,fr(unit))),spcounts);
end