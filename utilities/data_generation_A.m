function [W1,Wp,Clist_GT,threshold_GT, case_mtx, ctrl_mtx] = data_generation_A(cluster_size, case_num, ctrl_num, mu1, mu0, sigma1, sigma0, N, d_in_ratio, d_out_ratio)
% a simple data generator
% Chuan

% output: W1 ground truth
%         Wp permuted data
idx_temp = pick_idx(cluster_size, N);
% rng(0); 
pm = randperm(N);
N = N*(N-1)/2;cluster_size = length(idx_temp);
all_idx = 1:N;
% case_idx = randsample(node_num,cluster_size);
case_idx = idx_temp;
ctrl_idx = all_idx;
ctrl_idx(case_idx) = [];

% add non central dist. to the PT
d_out = floor((1-d_out_ratio) * length(ctrl_idx));
d_in = floor(d_in_ratio * length(case_idx));
d_out_idx = randsample(1:length(ctrl_idx),d_out,false);
d_out_orig = ctrl_idx(d_out_idx);
ctrl_idx(d_out_idx) = [];
d_in_idx = randsample(1:length(case_idx),d_in,false);
d_in_orig = case_idx(d_in_idx);
case_idx(d_in_idx) = [];


case_mtx = zeros(case_num, N);
case_mtx(:,case_idx) = randn(case_num, length(case_idx))*sigma1 + mu1;
case_mtx(:,d_out_orig) = randn(case_num, d_out)*sigma1 + mu1;
case_mtx(:,ctrl_idx) = randn(case_num, length(ctrl_idx))*sigma0 + mu0;
case_mtx(:,d_in_orig) = randn(case_num, d_in)*sigma0 + mu0;
ctrl_mtx = randn(ctrl_num, N)*sigma0 + mu0;

[~,p_vec]=ttest2(case_mtx,ctrl_mtx);
W1=squareform(-log(p_vec)); 
Wp = W1(pm,pm);
% return ground truth Clist
[~,Clist_GT] = sort(pm);
%% find a good cut
threshold_vec = [0.001 0.005 0.01 0.05 0.1];
f1score_vec = zeros(length(threshold_vec),1);
for i = 1:length(threshold_vec)
    target = ones(1,N);
    target(ctrl_idx)=0;
    output = p_vec <= threshold_vec(i);
    [~,cm,~,~] = confusion(target,output);
    f1score_vec(i) = cm(2,2)/(cm(2,2) + 0.5*(cm(1,2) + cm(2,1)));
end

[~,thresh_idx] = max(f1score_vec);
threshold_GT=threshold_vec(thresh_idx);
% sig_edge = sum(p_vec <= threshold_GT);

end