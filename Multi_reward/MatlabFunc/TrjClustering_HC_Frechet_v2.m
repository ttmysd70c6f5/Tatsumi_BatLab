function TrjCluster = TrjClustering_HC_Frechet_v2(Flights,varargin)
% INPUTS
% r_vector: position of the target bat. # samples x dimension
% f_vector: vector of being 1 if flying and 0 if not flying

%% Parameters
p = inputParser;
addRequired(p, 'Flights'); % Structure of extracted flights [1,bats]
addOptional(p, 'n_smp',10); % downsampling size
addOptional(p, 'cutoff',1.0); % cutoff value for clustering
addOptional(p, 'min_clus',3); % minimum number of trajectories to be considered as a cluster
addOptional(p, 'n_dim',3); % number of dimension to consider. 2 = only x and y, 3 = x, y, z
addOptional(p,'distance','Frechet'); % distance metric
addOptional(p, 'iprange',[0 1]); % interporation range: 0 = start, 1 = end of whole trajectory.

parse(p, Flights, varargin{:});

n_smp = p.Results.n_smp; % number of downsampled points
cutoff = p.Results.cutoff; % Cutoff value for clustering
min_clus = p.Results.min_clus; % Minimum amounts of trajectory in a cluster
n_dim = p.Results.n_dim; % dimension of trajectory (x,y,z)
distance_method = p.Results.distance; % distance metric
iprange = p.Results.iprange; % interporation range
if iprange(1) >= iprange(2)
    error('The 1st element should be smaller than the 2nd element.')
end
if iprange(1)<0 || iprange(2)>1
    error('iprange shuold be within [0 1].')
end


%% Main of function
for bb = 1:length(Flights)
    n_flight = Flights(bb).n_flight;
    % Original flights
    r_flight = cell(n_flight,1);
    for ff = 1:n_flight % landing idx
        r_flight{ff} = Flights(bb).r_flight{ff}(:,1:n_dim,bb); % coordinate of flight trajectory
    end

    % Downsampling of flights
    r_flight_ds = NaN(n_flight,n_smp*n_dim); % sampled position (# of flights, # of sampling * [x,y,z])
    for ff = 1:n_flight % landing idx
        len_flight = size(r_flight{ff},1); % length of flight
        % w_interp = 1 + fl_len_ts .* (0.5 + 0.5 .* (linspace(-n_smp/2,n_smp/2,n_smp)).^3 / (n_smp/2)^3); % weight for interpolation: sigmoid
        w_interp = linspace(len_flight*iprange(1)+1,len_flight*iprange(2),n_smp); % linear weight
        r_flight_ds(ff,:) = reshape(interp1(1:len_flight,r_flight{ff}, w_interp),[1,n_smp*n_dim]); % interporation
    end

    % Run hierarchical clustering based on Frechet distance
    X = r_flight_ds; % Input variable for clustering
    Y = [];
    pairs = nchoosek(1:n_flight,2);
    for n = 1:length(pairs)
        Y(n) = DiscreteFrechetDist(squeeze(X(pairs(n,1),:,:)),squeeze(X(pairs(n,2),:,:))); % Compute Frechet distance between trajectories
    end
    Z = linkage(Y,'single');
    clu_idx = cluster(Z,'Cutoff',cutoff,'Criterion','distance');
    
    [clu_sz, clu_list] = groupcounts(clu_idx);
    
    % Assign the small clusters to -1
    new_clu_idx = clu_idx;
    small_clus = clu_list(clu_sz < min_clus);
    new_clu_idx(ismember(new_clu_idx,small_clus)) = -1;
    [new_clu_sz,new_clu_list] = groupcounts(new_clu_idx);
    
    % Sort cluster according to the size
    [~,sort_idx] = sort(new_clu_sz,1,"descend");
    clu_sz_sorted = new_clu_sz(sort_idx);
    clu_list_sorted = new_clu_list(sort_idx);
    
    % Re-order cluster list so that the outliers come first
    if sum(clu_list_sorted==-1) > 0
        if clu_sz_sorted(clu_list_sorted==-1)>0 && find(clu_list_sorted==-1) ~= 1
            loc_outlier = find(clu_list_sorted==-1);
            clu_sz_sorted = clu_sz_sorted([loc_outlier,1:loc_outlier-1,loc_outlier+1:end],:);
            clu_list_sorted = clu_list_sorted([loc_outlier,1:loc_outlier-1,loc_outlier+1:end],:);
        end
    end
    
    % Reassign the cluster index according to the size of cluster
    if clu_sz_sorted(clu_list_sorted==-1)>0
        for i = 2:length(clu_list_sorted)
            clu = clu_list_sorted(i);
            new_clu = i-1 + 1000; % add a big number for this moment to avoid unexpected merging of clusters during updating
            new_clu_idx(new_clu_idx==clu) = new_clu;
            clu_list_sorted(i) = new_clu;
            clu_sz_sorted(i) = sum(new_clu_idx==new_clu);
        end
    else
        for i = 1:length(clu_list_sorted)
            clu = clu_list_sorted(i);
            new_clu = i + 1000;
            new_clu_idx(new_clu_idx==clu) = new_clu;
            clu_list_sorted(i) = new_clu;
            clu_sz_sorted(i) = sum(new_clu_idx==new_clu);
        end
    end
    new_clu_idx(new_clu_idx~=-1) = new_clu_idx(new_clu_idx~=-1)-1000;
    clu_list_sorted(clu_list_sorted~=-1) = clu_list_sorted(clu_list_sorted~=-1) - 1000;
    
    n_clu = length(clu_list_sorted); % number of final clusters
    
    % Calculate correlation and Frechet distance with mean of the cluster
    clu_corr = zeros(n_flight,1);
    clu_dist = zeros(n_flight,1);
    for ic = 1:n_clu
        clu = clu_list_sorted(ic);
        mean_path = squeeze(mean(r_flight_ds(new_clu_idx==clu,:,:),1));
        clu_locs = find(new_clu_idx==clu);
        for ic2 = 1:length(clu_locs)
            clu_loc = clu_locs(ic2);
            clu_corr(clu_loc) = mean(spdiags(corr(squeeze(r_flight_ds(clu_loc,:,:)),mean_path),0));
            [clu_dist(clu_loc),~] = DiscreteFrechetDist(squeeze(r_flight_ds(clu_loc,:,:)),mean_path);
        end
    end
    
    % mean of takeoff point for each cluster
    
    r_ds_mean = zeros(n_clu,n_smp*n_dim);
    for iclu = 1:n_clu % extract one cluster
        clu = clu_list_sorted(iclu);
        idx = new_clu_idx == clu; % index for the current cluster
        num_flights_clu = sum(idx); % number of flights contained in the current cluster
        r_flight_clu = r_flight_ds(idx,:,:); % downsampled positions for the extracted flights
        r_ds_mean(iclu,:) = reshape(squeeze(mean(r_flight_clu,1)),[1 n_smp*n_dim]);
    end
    r_toff_mean = r_ds_mean(:,1:n_smp:n_smp*n_dim); % averaged position of takeoffs for each cluster (# cluster x dimension)
    r_land_mean = r_ds_mean(:,n_smp:n_smp:n_smp*n_dim); % averaged position of landings for each cluster (# cluster x dimension)
    
      
    % Save to the output struct
    TrjCluster(bb).clu_idx      = new_clu_idx; % cluster index for each cluster
    TrjCluster(bb).clu_sz       = clu_sz_sorted; % count of flights in each clutser
    TrjCluster(bb).clu_list     = clu_list_sorted; % list of cluster
    TrjCluster(bb).n_clu        = n_clu; % number of final clusters
    TrjCluster(bb).r            = r_flight; % flight position without downsampling
    TrjCluster(bb).r_ds         = r_flight_ds; % downsampled flight
    TrjCluster(bb).clu_corr     = clu_corr; % correlation with mean of the cluster
    TrjCluster(bb).clu_dist     = clu_dist; % distance with mean of the cluster
    TrjCluster(bb).r_toff_mean  = r_toff_mean; % averaged position for takeoffs
    TrjCluster(bb).r_land_mean  = r_land_mean; % averaged position for landings
    TrjCluster(bb).r_ds_mean    = r_ds_mean; % averaged downsampled position
end

end
