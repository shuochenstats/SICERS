function [Cindx,CID,Clist,T]=kpartite(W,p0,fig, kmeans_iter)
%%%% This function is for parsimonious detector of Kpartite structure


%%% Inputs:
%%%%%  nlogp:   a 1 by n vector of t test p-values from the raw data. If
%%%%%           data contains k nodes, n=k*(k-1)/2

%%%%%  p0:      A threshold on the p-values. We only do the clustering on the
%%%%%           significant edges.
%%%%%  fig:     1 for plotting the quality-quantity curve and the
%%%%%           clustering result
%%%%%           0 for not plotting
%%%%%  kmeans_iter: number of iterations for total cluster K selection, the
%%%%%           precision of K and the time complexity will increase with kmeans_iter



%%% Outputs:
%%%%%  Cindx:    the cluster index of every non-isolated node
%%%%%  CID:     the cluster index of every cluster in a power descending
%%%%%           order. i.e. CID(1) will be the cluster index of the most
%%%%%           concentrated cluster
%%%%%  Clist:   the reordered node index, nodes in the same cluster are
%%%%%           permuted together in such way: [find(Cindx==CID(1))
%%%%%           find(Cindx==CID(2)) ... find(Cindx==CID(K))]

%%% Toy example: 
% % %     tau=0.3;
% % %     varval1=0.3; %varval1=0.1;
% % %     varval2=1; %varval2=0.1;
% % %     TP_edges=zeros(90,90);
% % %     TP_edges(1:10,11:20)=1;
% % %     TP_edges(11:20,1:10)=1;
% % %     for i=1:20
% % %         TP_edges(i,i)=0;
% % %     end
% % %     TP_edges=squareform(TP_edges);
% % %     inside_idx=find(TP_edges==1);
% % %     outside_idx=find(TP_edges==0);
% % %     casedata_cor=zeros(30,4005);
% % %     casedata_cor(:,inside_idx)=varval1*randn(30,length(inside_idx))+tau;
% % %     casedata_cor(:,outside_idx)=varval2*randn(30,length(outside_idx));
% % %     ctrldata_cor(:,inside_idx)=varval1*randn(30,length(inside_idx));
% % %     ctrldata_cor(:,outside_idx)=varval2*randn(30,length(outside_idx));
% % %     datacomb_cor=[casedata_cor;ctrldata_cor];
% % %     [h1,p1]=ttest(casedata_cor,ctrldata_cor); 
% % %     W1=squareform(-log(p1));
% % %     signet=W1(1:20,1:20);
% % %     ar=randperm(20);
% % %     signet=squareform(signet(ar,ar));
% % %     [Cindx,CID,Clist,T]=kpartite(signet,-log(0.1),1, 10);


%nlogp = -log(P)
%W is the locfdr score
W1=squareform(W);
W(W<p0)=0;%Threshold on the p-values
W=squareform(W);
 if fig==1   
figure;imagesc(W)
 end
z1=find(sum(W)>0);

W=W(z1,z1);

degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
L=D-W;%Laplacian matrix


[V,D]=eig(L);
%figure;plot(diag(D),'x')

diff   = eps;

lenW=length(find(W>0))/2;
sumW=sum(sum(W))/2;

%% Determine the number of clusters K
Mk=[];
Qual=[];
for m=1:kmeans_iter,
Prp_net=[];
Prp_net2=[];
for K=2:size(L,1)-1
   [U, ddd] = eigs(L,K);
C=kmeans(U,K);
indx=[]; %indx
A_net=[];% chi2 in the net
net_V=[];% size of each cluster
C_net=[];
for k=1:K
    %k=1
    indx=[indx;find(C==k)];
    net_V(k)=length(find(C==k));
    WC=W(find(C==k), find(C==k));
    W_off=W(find(C==k),setdiff(1:1:size(W,1),find(C==k)));
    %C_net(k)=length(find(W_off>0));
    C_net(k) = sum(W_off(find(W_off>0)))/2;
    A_net(k)=(net_V(k)*(net_V(k)-1))/2;
end
%Prp_net(K)=(sum(C_net)/sumW)*(sum(C_net)/(size(W,1)*(size(W,1)-1)/2-sum(A_net)));
%Prp_net(K)=(sum(C_net)/lenW)*(sum(C_net)/(size(W,1)*(size(W,1)-1)/2-sum(A_net))-(sumW-sum(C_net))/sum(A_net));
Prp_net2(K)= (sum(C_net)/(size(W,1)*(size(W,1)-1)/2-sum(A_net)));

end   


%K = find(Prp_net == max(Prp_net));

K = find(Prp_net2 == max(Prp_net2));
K=K(1);
Mk(m,:)=[K max(Prp_net2)];

Qual(:,m)=Prp_net2;
end
if fig==1
figure;plot(Qual,'x')
end
K=Mk(find(Mk(:,2)==max(Mk(:,2))),1);
K=K(1);
%% Find the cluster ID for each of the nodes
[U, ddd] = eigs(L,K-1);% zero eigenvalue makes this K-1
C=kmeans(U,K);

indx=[]; %indx
A_net=[];% chi2 in the net
net_V=[];% size of each cluster
A_net_off=[];
C_net=[];
C_net_w=[];
for k=1:K
    indx=[indx;find(C==k)];
    net_V(k)=length(find(C==k));
    WC=W(find(C==k), find(C==k));
    W_off=W(find(C==k),setdiff(1:1:size(W,1),find(C==k)));
    C_net(k)=length(find(W_off>0))/2;
    %C_net(k) = sum(W_off(find(W_off>0)))/2;
    C_net_w(k) = sum(W_off(find(W_off>0)))/2;
    A_net(k)=(net_V(k)*(net_V(k)-1))/2;
    A_net_off(k)=net_V(k)*(size(W,1)-net_V(k));
end
T=(sum(C_net)/lenW)*(sum(C_net)/(size(W,1)*(size(W,1)-1)/2-sum(A_net)));

   offdiagscore=(C_net_w).^2./(A_net_off);
   offdiagscore(isnan(offdiagscore))=0;
   [offdiagscore_sort,offdiagscore_sortID]=sort(offdiagscore,'descend');
   
inx_imporance=[];
for i=1:K,
   inx_imporance=[  inx_imporance; find(C==offdiagscore_sortID(i))];
end

Cindx = 1:size(W1,1);
Cindx(z1)=C;
Cindx(setdiff(1:size(W1,1),z1))=-1;
CID=offdiagscore_sortID;
Clist = z1(inx_imporance);
Clist = [Clist setdiff(1:size(W1,1),z1)];
if fig==1
figure; imagesc(W(Clist,Clist))
end