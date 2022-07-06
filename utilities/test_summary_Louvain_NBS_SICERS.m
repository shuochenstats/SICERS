function res_summary = test_summary_Louvain_NBS_SICERS(parameters, show, M,  bet_noise)
% This code summarizes 4 testing algorithms
%
% Input:
% parameters that generated the ground truth graph:
% mu0: mean value for the control group
% mu1: mean value for the case group
% sigma0: SD for the control group;
% sigma0: SD for the case group;
% cluster size = size of the only cluster
% ctrl_num: number cases in the of control group;
% case_num: number of cases in the case group;
% alpha parameter value under null hypothesis
% W1: Groud Truth
% Wp: Permuted (observed) data
% Clist_GT: Ground truth node list
% alpha: threshold for the -log(p0) table
% show: =0 or =1
% save: =0 or =1
% M = number of permutations in a permutation test

% Output:
% res_summary: structured data containing resulted recovered subgraphs,
% including the computed times, and all other variables that are needed,
% e.g. FDR, TPR, etc.

% Chuan Bi
% 01/31/2022

%% load parameters
mu0 = parameters.mu0; % default 0
mu1 = parameters.mu1;
sigma0 = parameters.sigma0; % default sigma0 = sigma1 = 1
sigma1 = parameters.sigma1;
cluster_size = parameters.cluster_size;
ctrl_num = parameters.ctrl_num; % fixed at 30
case_num = parameters.case_num; % fixed at 30
N = parameters.N;
d_in_ratio = parameters.d_in_ratio;
d_out_ratio = parameters.d_out_ratio;
if length(cluster_size) == 1
    [W1,Wp,Clist_GT,threshold_GT,Wt, Wt2] = data_generation_A(cluster_size, case_num, ctrl_num, mu1, mu0, sigma1, sigma0, N,d_in_ratio,d_out_ratio);
else
    [W1,Wp,Clist_GT,threshold_GT,Wt, Wt2] = data_generation_B(cluster_size, case_num, ctrl_num, mu1, mu0, sigma1, sigma0, N,d_in_ratio,d_out_ratio,bet_noise);
end
edge_density = @(C) sum(sum(C))/(size(C,1)*(size(C,1)-1));
node_density = @(C) sum(sum(C))/size(C,1);
custom_density = @(C) sum(sum(C))/(size(C,1));
cluster_nodes = {};
for i = 1:length(cluster_size)
    if i == 1
        cluster_nodes{i} = Clist_GT(1:cluster_size(1));
    else
        cluster_nodes{i} = Clist_GT(sum(cluster_size(1:i-1))+1:sum(cluster_size(1:i)));
    end
end
%% Dense subgraph detection greedy algorithm
tic
[W_greedy, Clist_greedy,CID_greedy] = greedy(Wp,  node_density, threshold_GT);
elapsed_greedy = toc;
%% Louvain Method
tic;
[W_Louvain, Clist_Louvain, CID_Louvain] = Louvain(Wp, node_density, threshold_GT);
elapsed_Louvain = toc;

%% NBS
tic;
threshold_NBS = -log(threshold_GT);
[W_NBS, Clist_NBS, CID_NBS] = NBS(Wp, threshold_NBS);
elapsed_NBS = toc;

%% SICERS
tic
[CID_SICERS,W_SICERS, Clist_SICERS]=SICERS_old(Wp,threshold_GT,10);
elapsed_SICERS = toc;

%% FDRBH
p_vec = exp(-squareform(Wp));
true_edge = zeros(N);
for i = 1:length(cluster_size)
    true_edge(cluster_nodes{i},cluster_nodes{i}) = 1;
end
true_edge = true_edge - diag(diag(true_edge));

false_edge = ones(N) - true_edge;
false_edge = false_edge - diag(diag(false_edge));
true_edge = squareform(true_edge);
false_edge = squareform(false_edge);

p1 = p_vec;
% [FDR] = mafdr(p1,'BHFDR',true);
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p1,0.2);
p1(p1>crit_p)=1;
    
W1_fdr=squareform(-log(p1));
% imagesc(W1_fdr);colormap jet;colorbar;
    
fdr_in_edge = squareform(W1_fdr)>0;
fdr_out_edge = 1-fdr_in_edge;

tpr_fdr = sum(true_edge.*fdr_in_edge)/sum(true_edge);
fdr_fdr =  sum(fdr_in_edge.*false_edge)/sum(fdr_in_edge);

%% summarize results

if show == 1
    figure;
    subplot(3,2,1)
    imagesc(W1);
    colormap jet;colorbar;
    title('Ground Truth','FontSize',30)
    
    subplot(3,2,2)
    imagesc(W1_fdr);
    colormap jet;colorbar;
    title('FDR-BH','FontSize',30)
    
    subplot(3,2,3)
    imagesc(W_NBS);
    colormap jet;colorbar;
    title('NBS','FontSize',30)
    
    subplot(3,2,4)
    imagesc(W_Louvain);colormap jet;colorbar;
    title('Louvain','FontSize',30)
    
    subplot(3,2,5)
    imagesc(W_SICERS);colormap jet;colorbar;
    title('SICERS','FontSize',30)
    
    subplot(3,2,6)
    imagesc(W_greedy);colormap jet;colorbar;
    title('greedy','FontSize',30)
    
    set(gcf,'position',[1001          46        1103        1302]);
    drawnow
end
%% save files
res_summary.parameters = parameters;
res_summary.Wp = Wp;
res_summary.W1 = W1;
res_summary.Clist_GT = Clist_GT;
res_summary.threshold_GT = threshold_GT;
res_summary.Wt = Wt;
res_summary.Wt2 = Wt2;
% res_summary.sig_edge = sig_edge;
% 
res_summary.greedy.W_greedy = W_greedy;
res_summary.greedy.Clist_greedy = Clist_greedy;
res_summary.greedy.elapsed_greedy = elapsed_greedy;
res_summary.greedy.CID_greedy = CID_greedy;


res_summary.Louvain.W_Louvain = W_Louvain;
res_summary.Louvain.Clist_Louvain = Clist_Louvain;
res_summary.Louvain.elapsed_Louvain = elapsed_Louvain;
res_summary.Louvain.CID_Louvain = CID_Louvain;

res_summary.NBS.W_NBS = W_NBS;
res_summary.NBS.Clist_NBS = Clist_NBS;
res_summary.NBS.elapsed_NBS = elapsed_NBS;
res_summary.NBS.CID_NBS = CID_NBS;
% 
res_summary.SICERS.W_SICERS = W_SICERS;
res_summary.SICERS.Clist_SICERS = Clist_SICERS;
res_summary.SICERS.elapsed_SICERS = elapsed_SICERS;
res_summary.SICERS.CID_SICERS = CID_SICERS;
% 
P_value_greedy = permutation_testA(res_summary, M, @greedy, node_density);
res_summary.greedy.P_value_Louvain = P_value_greedy;
P_value_Louvain = permutation_testA(res_summary, M, @Louvain, node_density);
res_summary.Louvain.P_value_Louvain = P_value_Louvain;
P_value_NBS = permutation_testA(res_summary, M, @NBS, node_density);
res_summary.NBS.P_value_NBS = P_value_NBS;
P_value_SICERS = permutation_testA(res_summary, M, @SICERS, node_density);
res_summary.SICERS.P_value_SICERS = P_value_SICERS;



%% determine false positives and FWER
% greedy
inference_greedy = eval_FD(cluster_nodes, Clist_greedy, CID_greedy, P_value_greedy);
res_summary.greedy.inference_greedy = inference_greedy;


% Louvain
inference_Louvain = eval_FD(cluster_nodes, Clist_Louvain, CID_Louvain, P_value_Louvain);
res_summary.Louvain.inference_Louvain = inference_Louvain;
% 
% NBS
inference_NBS = eval_FD(cluster_nodes, Clist_NBS, CID_NBS, P_value_NBS);
res_summary.NBS.inference_NBS = inference_NBS;
% 
% SICERS
inference_SICERS = eval_FD(cluster_nodes, Clist_SICERS, CID_SICERS, P_value_SICERS);
res_summary.SICERS.inference_SICERS = inference_SICERS;
% 
res_summary.FDRBH.power = tpr_fdr;
res_summary.FDRBH.FDR = fdr_fdr;

end