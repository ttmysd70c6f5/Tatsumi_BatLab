function [P,C] = comp_transmat(bpos,num_vars,n_tags,exc,type,random_test)
% This function computes transition matrix
% Input:
%     bpos: matrix labeled with spots for one bat
%     num_vars: number of states/spots
%     n_tags: number of bats
%     exc: values not to count as spots
%     type: P will contain transition matrix for spots when type='spot'
%         while matrix for transition between states will be generated when
%         type='state.'
%     random_test: true=randomize by circle shift
P = nan(num_vars,num_vars,n_tags); % transition matrix
C = zeros(num_vars,n_tags); % Count of states/spots

for i=1:n_tags
    bpos_shifted = bpos;
    % Extract the timings of landings and takeoffs
    [blanding,btakeoff] = findchanges(bpos_shifted,exc,1);
    if random_test % randomization test by shifting the positions of i-th bat
        IFI = diff(find(btakeoff(:,i)));
        shft_dist = [prctile(IFI,95),size(bpos_shifted,1)-prctile(IFI,95)]; % uniform distribution for circle shift
        shft = randi(int32(shft_dist),1);
        btakeoff(:,[1:i-1,i+1:end]) = circshift(btakeoff(:,[1:i-1,i+1:end]),shft,1);
        blanding(:,[1:i-1,i+1:end]) = circshift(blanding(:,[1:i-1,i+1:end]),shft,1);
        bpos_shifted(:,[1:i-1,i+1:end]) = circshift(bpos_shifted(:,[1:i-1,i+1:end]),shft,1);
    end
    
    % Identify the spots when bats takeoff and land
    bspot_land = bpos_shifted(blanding(:,i),:); 
    bspot_land = bspot_land(2:end,:); % remove the landing to initial spot
    bspot_toff = bpos_shifted(btakeoff(:,i),:);
    if ismember(bpos_shifted(end,i),exc) % if the i-th bat flied when task finished
        bspot_toff = bspot_toff(1:end-1,:); % remove the final takeoff
    end
    
    bspots = zeros(size(bspot_land,1),2,n_tags); % column 1 is takeoff spot, column 2 is landing spot for D2
    bspots(:,1,:) = bspot_toff; 
    bspots(:,2,:) = bspot_land;


    % Label the takeoff and landing spots
    bstate = zeros(size(bspots,1),2); % state transition vector for the i-th bat
    for j=1:length(bstate)
        for k=1:2 % takeoff and landing
            bspot_now = bspots(j,k,i); % current spot (takeoff/landing)
            switch type
                case 'state'   
                    % label current state
                    if ismember(bspot_now,1:100) && sum(bspots(j,k,:)==bspot_now)>=2
                        bstate(j,k) = 1; % social-home: home (1-3) and social (>= 2 animals)
                        C(1,i) = C(1,i)+1;
                    elseif ismember(bspot_now,1:100) && sum(bspots(j,k,:)==bspot_now)==1
                        bstate(j,k) = 2; % unsocial-home: home (1-3) and unsocial (< 2 animlas)
                        C(2,i) = C(2,i)+1;
                    elseif ismember(bspot_now,101:200) && sum(bspots(j,k,:)==bspot_now)>=2
                        bstate(j,k) = 3; % social-feeder: feeder (4-7) and social (>= 2 animals)
                        C(3,i) = C(3,i)+1;
                    elseif ismember(bspot_now,101:200) && sum(bspots(j,k,:)==bspot_now)==1
                        bstate(j,k) = 4; % unsocial-feeder: feeder (4-7) and unsocial (< 2 animals)
                        C(4,i) = C(4,i)+1;
                    else
                        error('Incorrect spot specified')
                    end
                case 'spot'   
                    % label current state
                    if ismember(bspot_now,1:100)
                        bstate(j,k) = bspot_now; % home (1:100)
                        C(bspot_now,i) = C(bspot_now,i)+1;
                    elseif ismember(bspot_now,101:200)
                        bstate(j,k) = bspot_now - 100; % home (101:200)
                        C(bspot_now-100,i) = C(bspot_now-100,i)+1;
                    else
                        error('Incorrect spot specified')
                    end
            end
        end
    end

    
    % In the case that zero probability is represented as NaN
    for j=1:size(bstate,1)
        if isnan(P(bstate(j,1),bstate(j,2),i))
            P(bstate(j,1),bstate(j,2),i) = 1;
        else
            P(bstate(j,1),bstate(j,2),i) = P(bstate(j,1),bstate(j,2),i) + 1;
        end
    end
    
end

end
