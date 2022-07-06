function [W1,Wp,Clist_GT,threshold_GT, case_mtx, ctrl_mtx] = data_generation_B(cluster_size, case_num, ctrl_num, mu1, mu0, sigma1, sigma0, N, d_in_ratio, d_out_ratio,bet_noise)
% a simple data generator
% This adds noisy signals on the inter-blocks region.
% Chuan

% output: W1 ground truth
%         Wp permuted data
%%
% initialization
[idx_temp, between_cluster_idx] = pick_idx(cluster_size, N);
case_mtx = zeros(case_num, N);
pm = randperm(N);
N_edge = N*(N-1)/2;
all_idx = 1:N_edge;
d_bet_ratio = 0.2;
%%
% assign mu2, sigma2 to the between
% rng(0); 
% case_idx = randsample(node_num,cluster_size);
case_idx = idx_temp;
ctrl_idx = all_idx;
ctrl_idx(case_idx) = [];

% add non central dist. to the PT

% assign 30% of the between cluster edges to N(mu1, sigma1)
d_bet = floor(d_bet_ratio * length(between_cluster_idx));
d_bet_idx = randsample(1:length(between_cluster_idx), d_bet, false);
between_cluster_idx =  between_cluster_idx(d_bet_idx);

d_out = floor((1-d_out_ratio) * length(ctrl_idx));
d_in = floor(d_in_ratio * length(case_idx));
d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
d_out_orig = ctrl_idx(d_out_idx);
ctrl_idx(d_out_idx) = [];
d_in_idx = randsample(1:length(case_idx),d_in,false);
d_in_orig = case_idx(d_in_idx);
case_idx(d_in_idx) = [];
% define patients' correlation matrices
case_mtx(:,case_idx) = randn(case_num, length(case_idx))*sigma1 + mu1;
case_mtx(:,d_out_orig) = randn(case_num, d_out)*sigma1 + mu1;
case_mtx(:,ctrl_idx) = randn(case_num, length(ctrl_idx))*sigma0 + mu0;
case_mtx(:,d_in_orig) = randn(case_num, d_in)*sigma0 + mu0;
if bet_noise == 1
    case_mtx(:,between_cluster_idx) = randn(case_num, length(between_cluster_idx))*sigma1 + mu1;
end
% define controlling subjects' correlation matrices
ctrl_mtx = randn(ctrl_num, N_edge)*sigma0 + mu0;

[~,p_vec]=ttest2(case_mtx,ctrl_mtx);
W1=squareform(-log(p_vec)); 
Wp = W1(pm,pm);
% return ground truth Clist
[~,Clist_GT] = sort(pm);
%% find a good cut
threshold_vec = [0.001 0.005 0.01 0.05 0.1];
f1score_vec = zeros(length(threshold_vec),1);
for i = 1:length(threshold_vec)
    target = ones(1,N_edge);
    target(ctrl_idx)=0;
    output = p_vec <= threshold_vec(i);
    [~,cm,~,~] = confusion(target,output);
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

[~,thresh_idx] = max(f1score_vec);
threshold_GT=threshold_vec(thresh_idx);
% sig_edge = sum(p_vec <= threshold_GT);

end