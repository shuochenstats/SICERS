%% SICERS network detection and testing
% Network object oriented statistics method is a new appraoch to detect and
% test hidden phenotype related connectivity patterns at the "network"
% level. 
% Notet that the subgraph may vary a little bit due to the random
% intialization of the kmeans++ algorithm. 
    %% D2 data analysis
    
    load('data_d2.mat')
    imagesc(WnTest);colorbar
    colormap jet
    title('Fig 1: D2 Heatmap of -log transformed p-value  matrix  ');
    snapnow;
    warning('off','all')
   [Test_Cindxn ,Test_CIDn ,Test_Clistn ]=SICERS_A(squareform(WnTest0),1,0,8); 
 

 

    %% Graph combinatorics based test 
    
    
    [signodeGEP,GEPstat,P_SICERS]=GEP_newstats_testonly(WnTest,WnTest0,Test_Cindxn,Test_CIDn,100);  
    
    P_SICERS
    
    
    imagesc(WnTest( Test_Clistn , Test_Clistn ));colorbar
    colormap jet
    title('Fig 2: Heatmap of -log transformed p-values with detected subnetworks');
    snapnow;