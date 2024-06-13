function [social_proximity, p_social_proximity, proximity_index] = compute_social_proximity(BhvData,n_ite)

r = BhvData.r;
n_bats = BhvData.n_tags;
bat_pairs = BhvData.bat_pairs; % all possible pairs
n_pair = size(bat_pairs,1); % number of pairs
% n_ite = 1000;

social_proximity = NaN(n_bats,n_bats,1+n_ite);
% p_social_proximity = NaN(n_bats);

for pair = 1:n_pair
% pair = 14; % current pair

    b1 = bat_pairs(pair,1);
    b2 = bat_pairs(pair,2);
    shft = randi(size(r,1),[n_ite,1]);
    shft = [0; shft];
    
    % uni-directional social proximity
    % dist_pair = [];
    % for ll = 1:length(Landings(b1).dist_oth_rest)
    %     dist_pair = [dist_pair;Landings(b1).dist_oth_rest{ll}(:,b2)];
    % end
    % social_proximity_b1_b2 = sum(dist_pair < 0.3) / size(dist_pair,1); % time fraction in which the distance between two bats was lower than 0.3m
    % 
    % dist_pair = [];
    % for ll = 1:length(Landings(b2).dist_oth_rest)
    %     dist_pair = [dist_pair;Landings(b2).dist_oth_rest{ll}(:,b1)];
    % end
    % social_proximity_b2_b1 = sum(dist_pair < 0.3) / size(dist_pair,1); % time fraction in which the distance between two bats was lower than 0.3m
    
    % bi-directional social proximity
    for ite = 1:1+n_ite
        dist_pair = vecnorm(r(:,:,b1)-circshift(r(:,:,b2),shft(ite)),2,2);
        social_proximity(b1,b2,ite) = sum(dist_pair < 0.3) / size(dist_pair,1); % time fraction in which the distance between two bats was lower than 0.3m
    end
end
social_proximity = social_proximity * 100;
p_social_proximity = sum(social_proximity(:,:,2:end) > social_proximity(:,:,1),3)/(size(social_proximity,3)-1);
proximity_index = social_proximity(:,:,1) - mean(social_proximity(:,:,2:end),3);

save(fullfile(pwd,'analysis','social_proximity.mat'),'social_proximity','p_social_proximity','proximity_index')