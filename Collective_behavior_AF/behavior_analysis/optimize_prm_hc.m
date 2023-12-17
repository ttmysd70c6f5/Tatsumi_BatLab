function [best_cutoff, num_good_clus] = optimize_prm_hc(r_interp, cutoffs, min_trj)

num_good_clus = zeros(length(cutoffs),1);
for i = 1:length(cutoffs)
    cutoff = cutoffs(i);
        % clustering
        hc_tree = linkage(r_interp,'single','euclidean'); % hierarchical cluster tree
        D = pdist(r_interp);
        leafOrder = optimalleaforder(hc_tree,D);
        hc_clus = cluster(hc_tree,'cutoff',cutoff,'Criterion','inconsistent');
        [clus_cts,clus_id] = groupcounts(hc_clus);
        clus_good = clus_id(clus_cts > min_trj);

        num_good_clus(i) = length(unique(clus_good));
end

[~,best_cutoff_id] = max(num_good_clus); % Extract the best cutoff value
best_cutoff = cutoffs(best_cutoff_id);
