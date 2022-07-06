function [W_DSD_greedy, Clist,Node_Seq,  removing_node] = OQC_greedyA(Wp_DSD,func)
% a simple implementation of the greedy algorithm from the Denser than
% dense paper
% note that this only extract ONE dense subgraphs
% Chuan

% Wp_DSD(Wp_DSD<threshold_DSD) = 0;
N = size(Wp_DSD,1);
Recording_Matrix = [];
Recording_Clist = 1:N;
Wp_temp = Wp_DSD;
for i = N:-1:1
    idxlist_temp = 1:length(Wp_temp);
    % work on the temp matrix
    [~,idx_min_temp] = min(sum(Wp_temp));
    % find corresponding index for idx_min_temp in Recording_Clist
    idxlist_temp(idx_min_temp) = [];
    
    score_temp = func(Wp_temp(idxlist_temp,idxlist_temp));
    % record index and score
    Recording_Matrix = [Recording_Matrix; [Recording_Clist(idx_min_temp),score_temp] ];
    Recording_Clist(idx_min_temp) = [];
    % update matrix
    Wp_temp = Wp_temp(idxlist_temp,idxlist_temp);
end
[~,max_idx] = max(Recording_Matrix(:,2));
removing_node = Recording_Matrix(1:max_idx,1);
Node_Seq = Recording_Matrix(end:-1:max_idx+1,1);
Clist = [Node_Seq;removing_node];
W_DSD_greedy = Wp_DSD([Node_Seq;removing_node],[Node_Seq;removing_node]);

end