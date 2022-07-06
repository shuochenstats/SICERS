function [W_Louvain, Clist_Louvain,CID] = Louvain(Wp, func, threshold_Louvain)
% Louvain method
% added a threshold to be the mean value of the Wp to determine whether to
% consider a community detected by Louvain or not.


Wp_Louvain = Wp;
Wp_Louvain(Wp_Louvain<-log(threshold_Louvain)) = 0;
mean_Wp = mean(mean(Wp_Louvain));


[Cindex,~]=community_louvain(Wp_Louvain);
coms = sort(unique(Cindex));
temp_Louvain = {};
res_Louvain = [];
blk_Size = [];
non_informative_idx = [];
I = 1;
for i = 1:length(unique(Cindex))
    idx = find(Cindex == coms(i));
    % tell if a community has greater mean values than the threshold
    if mean(mean(Wp_Louvain(idx,idx))) > mean_Wp
        blk_Size = [blk_Size;length(idx)];
        temp_Louvain{I} = idx;
        res_Louvain = [res_Louvain; func(Wp_Louvain(idx,idx))];
        I = I + 1;
    else
        % get noninformative portion
        non_informative_idx = [non_informative_idx; idx];
    end
end
[~,sorted_idx] = sort(res_Louvain,'descend');
Clist_Louvain = [];
for i = 1:length(sorted_idx)
    Clist_Louvain = [Clist_Louvain;temp_Louvain{sorted_idx(i)}];
end
Clist_Louvain = [Clist_Louvain; non_informative_idx];
CID = blk_Size(sorted_idx);
W_Louvain = Wp(Clist_Louvain,Clist_Louvain);
end