function [idx_cluster, between_cluster_idx] = pick_idx(cluster_size, N)
% given a vector of cluster sizes, return the corresponding indices the
% vectorized edge list.


Z = zeros(N);

for i = 1:length(cluster_size)
    if i == 1
        Z(1:cluster_size(i),1:cluster_size(i)) = 1;
    else
        Z(sum(cluster_size(1:i-1))+1:sum(cluster_size(1:i)), sum(cluster_size(1:i-1))+1:sum(cluster_size(1:i))) = 1;
    end
end
Z_temp = Z(1:sum(cluster_size), 1:sum(cluster_size));
bet_blk_mask = Z_temp ~= 1; bet_blk_mask = bet_blk_mask * 2;
for i = 1:N
    Z(i,i) = 0;
end
idx_cluster = find(squareform(Z));

% find between cluster
Z = zeros(N);
Z(1:sum(cluster_size), 1:sum(cluster_size)) = Z(1:sum(cluster_size), 1:sum(cluster_size)) + bet_blk_mask;
between_cluster_idx = find(squareform(Z));
end