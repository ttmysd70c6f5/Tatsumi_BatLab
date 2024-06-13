function plot_reward_distance(d_vectors,rate_vector)
% d_vectors: col1 = rw 1, col2 = rw2
% rate_vector: firing rate

c_vector = ['r','b'];
hold on
for rw = 1:2 % reward type
    [d_vector, sort_idx] = sort(d_vectors(:,rw));
    rate_vector = rate_vector(sort_idx);
    
    edges = 0:0.2:7;
    bin_value = mean([edges(1:end-1);edges(2:end)],1)';
    [~,~,bin_vector] = histcounts(d_vector,edges);
    d_vector_binned = bin_value(unique(bin_vector));
    rate_vector_binned = groupsummary(rate_vector,bin_vector,'mean');
    sem_vector_binned = groupsummary(rate_vector,bin_vector,'std') ./ sqrt(groupsummary(rate_vector,bin_vector,'nnz'));
    
    boundedline(d_vector_binned,rate_vector_binned,sem_vector_binned,c_vector(rw),'alpha','transparency', 0.1)
end
ylabel('Firing rate (Hz)')
xlabel('Distance to feeder (m)')
set(gca,'FontSize',16)