
%% SICERS network detection and testing
% Network object oriented statistics method is a new appraoch to detect and
% test hidden phenotype related connectivity patterns at the "network"
% level. 
% Notet that the subgraph may vary a little bit due to the random
% intialization of the kmeans++ algorithm. 
    %% D1 data analysis
    
    load('data_d1.mat')
    imagesc(WnTr);colorbar
    colormap jet
    title('Fig 1: Heatmap of -log transformed p-values of 90*90 matrix');
    snapnow;
    warning('off','all')
    [Tr_Cindxn ,Tr_CIDn ,Tr_Clistn ]=SICERS_A(squareform(WnTr0),1,0,8);
     
    
    %% Graph Combinatorics based test
    [signodeGEP,GEPstat,P_SICERS]=GEP_newstats_testonly(WnTr,WnTr0,Tr_Cindxn,Tr_CIDn,100);  
    
    P_SICERS
    
    
    imagesc(WnTr( Tr_Clistn , Tr_Clistn ));colorbar
    colormap jet
    title('Fig 2: Heatmap of -log transformed p-values of 90*90 matrix with detected subnetworks');
    snapnow;