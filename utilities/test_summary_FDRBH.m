function res_summary = test_summary_FDRBH(parameters, show, M,  bet_noise)
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
p_vec = exp(-squareform(Wp));

%%
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
% fpr_fdr = sum(false_edge.*fdr_in_edge)/sum(false_edge);
% tnr_fdr = sum(false_edge.*fdr_out_edge)/sum(false_edge);
% fnr_fdr = sum(true_edge.*fdr_out_edge)/sum(true_edge);
fdr_fdr =  sum(fdr_in_edge.*false_edge)/sum(fdr_in_edge);

%% save files
res_summary.parameters = parameters;
res_summary.Wp = Wp;
res_summary.W1 = W1;
res_summary.Clist_GT = Clist_GT;
res_summary.threshold_GT = threshold_GT;
res_summary.Wt = Wt;
res_summary.Wt2 = Wt2;
res_summary.FDRBH.power = tpr_fdr;
res_summary.FDRBH.FDR = fdr_fdr;
end