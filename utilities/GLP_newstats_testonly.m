
function [signode,GLPstat,P_value]=GLP_newstats_testonly(Wt,Wt2,Cindx,CID,M)
% M is the number of label permutation


[h1,p1]=ttest2(Wt',Wt2');
Wt_c=[Wt Wt2];
W1=squareform(-log(p1));
Cindx0=Cindx;
CID0=CID;
K_0=size(CID,2); %number of clusters
allinx=1:size(W1,1);

T_k0=zeros(K_0,1);
k_0=zeros(K_0,1);% number of nodes in kth cluster
for c=1:K_0
    idx_c=find(Cindx==c);
    nidx_c=size(idx_c,2);
    k_0(c)=nidx_c;
    nedges=nidx_c*(nidx_c-1)/2;
    G_k=W1(idx_c,idx_c);
    if(nidx_c<=4)
        T_k0(c)=0;
    else
%         T_k0(c)=sqrt(sum(sum(G_k))/2/nedges);
%         z=sum(sum(G_k))/2/nedges;
%         T_k0(c)=-log((z*exp(1-z))^(nedges));
        xbar=sum(squareform(G_k))/nedges;
        T_k0(c)=nedges*(xbar-1-log(xbar));
    end
    close all;
end


T_max = zeros(M,1);
for m=1:M
%     ar=randperm(90);
%     Z90_1=Z90_1(ar,ar);
%     Z90_2=Z90_2(ar,ar);
%     aridx=arrayfun(@(x) find(x==ar),1:1:90);
%     aridx1_20=aridx(1:20);%truly connected nodes index in permutated ar


    ar=randperm(60);
    
    Wt_c=Wt_c(:,ar);
    Wt1_c=Wt_c(:,1:30);
    Wt2_c=Wt_c(:,31:60);
    [h1,p1]=ttest(Wt1_c',Wt2_c');
    
    W1=squareform(-log(p1));

    nlogp=squareform(W1);
    [Cindx,CID,Clist]=SICERS_A(nlogp,0.08,0,3);
    K=size(CID,2); %number of clusters
    G=squareform(nlogp);
    T_k=zeros(K,1);
    for c=1:K
        idx_c=find(Cindx==c);
        nidx_c=size(idx_c,2);
        nedges=nidx_c*(nidx_c-1)/2;
        G_k=G(idx_c,idx_c);
        if(nidx_c<=4)
            T_k(c)=0;
        else
%             z=sum(sum(G_k))/2/nedges;
%             T_k(c)=-log((z*exp(1-z))^(nedges));
            % The new statistic is max{|E|*(xbar-1-log(xbar))}
            % where xbar=-2\sum(log(p_ij))/2|E|
            xbar=sum(squareform(G_k))/nedges;
            T_k(c)=nedges*(xbar-1-log(xbar));
        end
        close all;
    end
    T_max(m)=max(T_k);
    
end

Tmax_5prct = prctile(T_max,95);
[diff_value,diff_ind]=min(abs(sort(T_max)-max(T_k0)));
P_value=1-diff_ind/M;


%% account for multiple significant clusters
signode={};
SCinx = find(T_k0>Tmax_5prct);
noc = size(SCinx,1);
GLPstat=zeros(noc,1);
for j=1:noc
    signode{j}=allinx(ismember(Cindx0,SCinx(j)));
    GLPstat(j) = T_k0(SCinx(j));
end
end