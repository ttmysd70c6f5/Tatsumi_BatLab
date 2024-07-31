function SocialFlights = AdjustACdistribution(SocialFlights,observe_bat,bin,seed)

n_bat = length(SocialFlights);

for flying_bat = [1:observe_bat-1,observe_bat+1:n_bat]
    % bname = SocialFlights(flying_bat).bname;
    clu_list = SocialFlights(flying_bat).clu_list;
    n_clu = length(clu_list);
    ac_adjusted = cell(n_clu,2); % social, asocial
    for clu = 1:n_clu
        idx_clu_now = SocialFlights(flying_bat).clu_idx == clu_list(clu);
        idx_social = SocialFlights(flying_bat).social & ~SocialFlights(flying_bat).obs_at_fd & SocialFlights(flying_bat).stationary & SocialFlights(flying_bat).is_include;
        idx_asocial = ~SocialFlights(flying_bat).social & ~SocialFlights(flying_bat).obs_at_fd & SocialFlights(flying_bat).stationary & SocialFlights(flying_bat).is_include;
        
        ac_social_idx = find(idx_clu_now & idx_social);
        ac_social = SocialFlights(flying_bat).a_abs_observe(ac_social_idx); % absolute acceleration of observor bat
        
        ac_asocial_idx = find(idx_clu_now & idx_asocial);
        ac_asocial = SocialFlights(flying_bat).a_abs_observe(ac_asocial_idx); % absolute acceleration of observor bat

        % n_flight_social = length(ac_social);
        % n_flight_asocial = length(ac_asocial);

        if ~isempty(ac_social) && ~isempty(ac_asocial)
            [ac_social_adjusted,ac_asocial_adjusted,ac_social_adjusted_idx,ac_asocial_adjusted_idx] = adjust_distribution(ac_social,ac_asocial,bin,seed);
            ac_adjusted{clu,1} = ac_social_adjusted;
            ac_adjusted{clu,2} = ac_asocial_adjusted;
        end
        
    end

    SocialFlights(flying_bat).a_adjusted = ac_adjusted;
end