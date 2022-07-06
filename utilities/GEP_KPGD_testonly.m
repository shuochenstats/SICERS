
function [signode,GLPstat,P_value]=GEP_KPGD_testonly(Wt,Wt2,Cindx,CID,p0,M)
% M is the number of label permutation


[h1,p1]=ttest2(Wt',Wt2');
Wt_c=[Wt Wt2];
W1=squareform(-log(p1));
Cindx0=Cindx;
CID0=CID;
K_0=size(CID,2); %number of clusters

T_k0=zeros(K_0,1);
k_0=zeros(K_0,1);% number of nodes in kth cluster
sum_diag=0;
nedges_diag=0;
for c=1:K_0
    idx_c=find(Cindx==c);
    nidx_c=size(idx_c,2);
    k_0(c)=nidx_c;
    nedges=nidx_c*(nidx_c-1)/2;
    G_k=W1(idx_c,idx_c);
    sum_diag=sum_diag+sum(sum(squareform(G_k)));
    nedges_diag=nedges_diag+nedges;
end

nedges_offdiag=length(p1)-nedges_diag;
sum_offdiag=sum(-log(p1))-sum_diag;
Tk0=sum_offdiag/nedges_offdiag-sum_diag/nedges_diag;



T_max = zeros(M,1);
for m=1:M
%     ar=randperm(90);
%     Z90_1=Z90_1(ar,ar);
%     Z90_2=Z90_2(ar,ar);
%     aridx=arrayfun(@(x) find(x==ar),1:1:90);
%     aridx1_20=aridx(1:20);%truly connected nodes index in permutated ar

    ar=randperm(size(Wt,1));        %permutation of edges
    
    Wt_c=Wt_c(ar,:);
    Wt1_c=Wt_c(:,1:30);
    Wt2_c=Wt_c(:,31:60);
    [h1,p1]=ttest(Wt1_c',Wt2_c');
    
    W1=squareform(-log(p1));

    nlogp=squareform(W1);
    [Cindx,CID,Clist,T]=kpartite(nlogp,p0,0,5);
    K=size(CID,2); %number of clusters
    G=squareform(nlogp);
    T_k=zeros(K,1);
    sum_diag=0;
    nedges_diag=0;
    for c=1:K
        idx_c=find(Cindx==c);
        nidx_c=size(idx_c,2);
        nedges=nidx_c*(nidx_c-1)/2;
        G_k=G(idx_c,idx_c);
        sum_diag=sum_diag+sum(sum(squareform(G_k)));
        nedges_diag=nedges_diag+nedges;
    end
    nedges_offdiag=length(p1)-nedges_diag;
    sum_offdiag=sum(-log(p1))-sum_diag;
    T_max(m)=sum_offdiag/nedges_offdiag-sum_diag/nedges_diag;
    
end

Tmax_5prct = prctile(T_max,95);
[diff_value,diff_ind]=min(abs(sort(T_max)-Tk0));
P_value=1-diff_ind/M;

%% account for multiple significant clusters
if Tk0>Tmax_5prct
    signode=1;
else
    signode=0;
end
GLPstat=Tk0;
end