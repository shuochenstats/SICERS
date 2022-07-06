function P_value = permutation_testA(res_summary, M, funcA, funcB)


% edge_density = @(C) sum(sum(C))/(size(C,1)*(size(C,1)-1));
% node_density = @(C) sum(sum(C))/size(C,1);

% switch between functions
Wt = res_summary.Wt; Wt2 = res_summary.Wt2;
% sig_edge= res_summary.sig_edge;
casenum = size(Wt,1);
ctrlnum = size(Wt2,1);
% Wp=res_summary.Wp;
Wt_orig = [res_summary.Wt;res_summary.Wt2];
threshold_GT = res_summary.threshold_GT;
Wp_orig = res_summary.Wp;
% determine which function to work with
cur_func = func2str(funcA);
cur_struct = res_summary.(cur_func);
CID0 = cur_struct.(['CID_',cur_func]);
Clist0 = cur_struct.(['Clist_',cur_func]);
r_vec = [0.1 0.05 0.01 0.005 0.001];
if strcmp(cur_func,'NBS')
    T_orig =  custom_statistic_NBS(Wp_orig, Clist0, CID0,-log(threshold_GT));
else
    T_orig =  custom_statistic(Wp_orig, Clist0, CID0, r_vec);
end


T_vec = zeros(M,1);
for m=1:M
    
    Wt_c=Wt_orig(randperm(casenum+ctrlnum),:);
    Wt1_c=Wt_c(1:casenum,:);
    Wt2_c=Wt_c(casenum+1:end,:);
    [~,ptemp]=ttest2(Wt1_c,Wt2_c);
    
    Wp=squareform(-log(ptemp));
    
    % select sig_edge edges such that the comparison is fair
%     [~,klg_idx] = maxk(-ptemp,sig_edge);
%     maxp = max(ptemp(klg_idx));
    if strcmp(cur_func,'NBS')
        [W_NBS, Clist, CID] = NBS(Wp, -log(threshold_GT));
        [~,T_vec(m)] = custom_statistic_NBS(Wp, Clist, CID, -log(threshold_GT));
    elseif strcmp(cur_func,'greedy')
        [W_greedy, Clist, CID] = greedy(Wp, funcB, threshold_GT);
        [~,T_vec(m)] = custom_statistic(Wp, Clist, CID, r_vec);
%         T_vec(m) = custom_statisticB(W_greedy, Clist, CID);
    elseif strcmp(cur_func,'Pard')
        [CID,W_Pard,Clist]=Pard_old(Wp,threshold_GT,5);
        [~,T_vec(m)] = custom_statistic(Wp, Clist, CID, r_vec);
%         T_vec(m) = custom_statisticB(W_Pard, Clist, CID);
    elseif strcmp(cur_func,'Louvain')
        [W_Louvain, Clist, CID] = Louvain(Wp, funcB, threshold_GT);
        [~,T_vec(m)] = custom_statistic(Wp, Clist, CID, r_vec);
%         T_vec(m) = custom_statisticB(W_Louvain, Clist, CID);
    end
    
    
    
end

P_value = zeros(length(T_orig),1);
for i = 1:length(T_orig)
    if CID0(i) <= 2
        P_value(i) = 1;
    else
        P_value(i) = sum((T_vec - T_orig(i)) >0 )/M;
    end
end




end