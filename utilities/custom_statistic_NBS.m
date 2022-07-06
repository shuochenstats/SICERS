function [T0_vec,score_max] = custom_statistic_NBS(Wp, Clist, CID, threshhold_NBS)
% thresholding
Wp_NBS = Wp;
Wp_NBS(Wp_NBS<threshhold_NBS) = 0;
Wp_NBS(Wp_NBS>=threshhold_NBS) = 1;
% get number of clusters
cluster_len = length(CID);
% define test statistic for each cluster
T0_vec = zeros(cluster_len,1);
% loop each cluster and calculate the statistic for each cluster
for i = 1:cluster_len
    if i == 1
        W_cluster = Wp_NBS(Clist(1:CID(1)),Clist(1:CID(1)));
    else
        W_cluster = Wp_NBS(Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)),...
            Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i)));
    end
    % calculate statistic for NBS, i.e. find the number of significant
    % edges
    if size(W_cluster,1) <= 5
        cluster_vec = 0;
    else
        cluster_vec = squareform(W_cluster);
    end
    T0_vec(i) = sum(cluster_vec);

end
% get the maximum of the scores across all clusters
score_max = max(T0_vec);

end