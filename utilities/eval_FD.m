function ret_struct = eval_FD(cluster_nodes, Clist, CID, P_value)
% evaluate False Discovery given the gound truth, denoted as cluster_nodes,
% and recovered clusteres, characterized as Clist and CID
% The P_value is a vector returned by the permutation test, to indicate
% whether a recovered cluster is significant or not.

% if there are more than 1 clusters in gound truth and returned more than 1
% clusters, we compare each returned cluster with each ground truth
% cluster, and return the one with the largest Jaccard index.

%% Network-wise inference
% the confusion matrix may look like
%                        P < 0.05       |   P >= 0.05
% ---------------------------------------------
%  J>= 50%|      TP                |    FN
% ---------------------------------------------
%   J < 50%|       FP               |    TN
% ---------------------------------------------

%% Edge-level inference
% Among the TP tests, we calculate the edge-wise inference of the subgraphs
% the confusion matrix may look like
%                             (i,j) in Gc hat  |  (i,j) not in Gc hat
% ----------------------------------------------------
% (i,j) in Gc        |      TP                |    FN
% -----------------------------------------------------
% (i,j) not in Gc |       FP               |    TN
% ----------------------------------------------------

% Chuan

%% Network inference first
TP_ntwk = 0;
FN_ntwk = 0;
FP_ntwk = 0;
TN_ntwk = 0;

 % for edge-wise inference
nodeLen = length(Clist);
I = 1;
for i = 1:length(P_value)
    Zero_GT = zeros(nodeLen);
    Zero_est = zeros(nodeLen);

    % load cluster node list
    if i == 1
        cur_list = Clist(1:CID(1));
    else
        cur_list = Clist(sum(CID(1:i-1))+1:sum(CID(1:i-1))+CID(i));
    end
    % calculate Jaccard Index
    JD_vec = zeros(length(cluster_nodes),1);
    for j = 1:length(cluster_nodes)
        cur_cluster = cluster_nodes{j};
        JD_vec(j) = length(intersect(cur_cluster,cur_list))/length(union(cur_cluster,cur_list));
    end
    [JD, imax_JD] = max(JD_vec);
    % we have found a cluster that matches the returned cluster closely,
    % call it cNodeTemp
    cNodeTemp = cluster_nodes{imax_JD};
    % first tell if this is significant
    if P_value(i) < 0.05
        if JD >= 0.5
 %% do edge-wise inference here
            Zero_GT(cNodeTemp,cNodeTemp) = 1; Zero_GT = Zero_GT - diag(diag(Zero_GT)); edge_GT = squareform(Zero_GT);
            Zero_est(cur_list,cur_list) = 1; Zero_est = Zero_est - diag(diag(Zero_est)); edge_est = squareform(Zero_est);
            % 2-by-2 confusion matrix, the structure looks like
    % ----------------------------------------------------------------
    %                                                     Pred
    % ----------------------------------------------------------------
    %                     |                    0                  1
    % -------------------------------------------------------------
    %        True    |    0 
    %                     |    1
            [c,cm,ind,per] = confusion(edge_GT,edge_est);
            ret_struct.(['cluster_',num2str(I)]).TPR = cm(2,2)/(cm(2,1) + cm(2,2));
            ret_struct.(['cluster_',num2str(I)]).FDR = cm(1,2)/(cm(1,2) + cm(2,2));
            ret_struct.(['cluster_',num2str(I)]).FPR = cm(1,2)/(cm(1,2) + cm(1,1));
            ret_struct.(['cluster_',num2str(I)]).FNR = cm(2,1)/(cm(2,1) + cm(2,2));

            % update 
            TP_ntwk = TP_ntwk + 1;
             ret_struct.(['cluster_',num2str(I)]).node_list= cur_list;
            I = I+1;

        else
            FP_ntwk = FP_ntwk + 1;
        end
    else
        if JD >= 0.5
            FN_ntwk = FN_ntwk + 1;
        else
            TN_ntwk = TN_ntwk + 1;
        end
    end
end

ret_struct.TP_ntwk = TP_ntwk;
ret_struct.TN_ntwk = TN_ntwk;
ret_struct.FP_ntwk = FP_ntwk;
ret_struct.FN_ntwk = FN_ntwk;









end
