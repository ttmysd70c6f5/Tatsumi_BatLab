function sp_mat = EventAlignedSpikes(SocialFlights,spikeData,bbb,clu,w_step,w_sz,t_margin,fr_th)

pts_margin = ceil(w_sz/w_step); % marginal points (remove later)

% time bin
t_bin = t_margin(1)-pts_margin*w_step:w_step:t_margin(2)+pts_margin*w_step;
t_edge = [t_bin(1)-w_step/2, t_bin+w_step/2];

% Spike data
target_units =  find(spikeData.good_units.fr >= fr_th); % Exclusion criteria
n_unit = length(target_units);

clu_idx = SocialFlights(bbb).clu_idx;

% Check exclusion criteria
is_social = SocialFlights(bbb).social(clu_idx==clu);
n1 = sum(is_social);
n2 = sum(~is_social);
sp_mat = NaN(n1+n2,length(t_bin),n_unit); % event aligned spikes
if n1 > 5 && n2 > 5 % minimum flights criteria

    % Event-aligned spike activity
    t_flight = SocialFlights(bbb).t_sec(clu_idx==clu);
    
    for unit = 1:n_unit
        sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{target_units(unit)} / 1e6; % spike times of the current unit in sec
    
        for ff = 1:n1+n2
            sp_mat(ff,:,unit) = histcounts(sptimes_unit_sec,t_flight(ff)+t_edge);
        end
    end
end