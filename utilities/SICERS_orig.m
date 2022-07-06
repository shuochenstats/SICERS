function [Cindx,CID,Clist]=SICERS_orig(nlogp,p0,fig, kmeans_iter)
%%%% This function is for parsimonious detector


%%% Inputs:
%%%%%  nlogp:   a 1 by n vector of t test p-values from the raw data. If
%%%%%           data contains k nodes, n=k*(k-1)/2

%%%%%  p0:      A threshold on the p-values. We only do the clustering on the
%%%%%           significant edges.

%%% Outputs:
%%%%%  Cindx:    the cluster index of every non-isolated node
%%%%%  CID:     the cluster index of every cluster in a power descending
%%%%%           order. i.e. CID(1) will be the cluster index of the most
%%%%%           concentrated cluster
%%%%%  Clist:   the reordered node index, nodes in the same cluster are
%%%%%           permuted together in such way: [find(Cindx==CID(1))
%%%%%           find(Cindx==CID(2)) ... find(Cindx==CID(K))]

 
%% Preprocessing of data

%nlogp = -log(P)
W=squareform(nlogp);
W1=W;
W(W1<-log(p0))=0;%Threshold on the p-values
    
%figure;imagesc(W)
z1=find(sum(W)>0); %Exclude the isolated nodes
% 
% if(isempty(z1))
%     % if after the screening the matrix is all zero
%     Cindx = ones(1,size(W1,1));
%     CID=1;
%     Clist = 1:size(W1,1);
% else

W=W(z1,z1);


degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
L=D-W;%Laplacian matrix


[V,D]=eig(L);
%figure;plot(diag(D),'x')

diff   = eps;
 
lenW=length(find(W>0))/2;

%% Determine the number of clusters K
Mk=[];
Qual=[];
for m=1:kmeans_iter,
    Prp_net=[];
    for K=1:size(L,1)
        try
            [U, ddd] = eigs(L,K, diff);
        catch
            [U, ddd] = eigs(L,K);
        end
        
        C=kmeans(U,K,'Replicates',kmeans_iter*2);
       
        indx=[]; %indx
        A_net=[];% chi2 in the net
        net_V=[];% size of each cluster
        C_net=[];
        for k=1:K
            indx=[indx;find(C==k)];
            net_V(k)=length(find(C==k));
            WC=W(find(C==k), find(C==k));
            %     C_net(k)=length(find(WC>0))/2;
            C_net(k) = sum(WC(find(WC>0)))/2;
            A_net(k)=(net_V(k)*(net_V(k)-1))/2;
        end
        Prp_net(K)=(sum(C_net)/lenW)*(sum(C_net)/sum(A_net));
    end
    
    
    
    K = find(Prp_net == max(Prp_net));
    K = K(1);%In case several k's give the same Prp_net value
    Mk(m,:)=[K max(Prp_net)];
    
    Qual(:,m)=Prp_net;
end
if fig==1
    figure;plot(Qual,'x')
end
K=Mk(find(Mk(:,2)==max(Mk(:,2))),1);
K=K(1);

%% Find the cluster ID for each of the nodes
try
    [U, ddd] = eigs(L,K, diff);
catch
   [U, ddd] = eigs(L,K);
end
C=kmeans(U,K);
indx=[]; 
A_net=[];
net_V=[];
C_net=[];
for k=1:K
    indx=[indx;find(C==k)];
    net_V(k)=length(find(C==k));
    WC=W(find(C==k), find(C==k));
%     C_net(k)=length(find(WC>0))/2;
    C_net(k) = sum(WC(find(WC>0)))/2;
    A_net(k)=(net_V(k)*(net_V(k)-1))/2;
end


diagscore=(C_net).^2./(A_net)/lenW;
diagscore(isnan(diagscore))=0;
[diagscore_sort,diagscore_sortID]=sort(diagscore,'descend');
   
inx_imporance=[];
for i=1:K,
   inx_imporance=[  inx_imporance; find(C==diagscore_sortID(i))];
end

Cindx = 1:size(W1,1);
Cindx(z1)=C;
Cindx(setdiff(1:size(W1,1),z1))=-1;
CID=diagscore_sortID;
Clist = z1(inx_imporance);
Clist = [Clist setdiff(1:size(W1,1),z1)];

end
%