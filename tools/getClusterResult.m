% function [mean_ACC, mean_NMI,mean_Purity,mean_Fscore,mean_Preci,mean_Recall,mean_ARI] = getClusterResult(Q,numClust,truthF,T,Dataname)
function [results,results1,results2,input] = getClusterResult(Q,Q1,Q2,numClust,truthF,T,Dataname,input)
if strcmp(Dataname,'handwritten-5view')
    for iter_c=1:T %
        C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        result_LatLRR = ClusteringMeasure(truthF, C);
        %         C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        %         result_LatLRR = EvaluationMetrics(truthF, C);  % res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];
        %c=1;
        acc(iter_c)    = result_LatLRR(1)*100;
        nmi(iter_c) = result_LatLRR(2)*100;
        purity(iter_c)= result_LatLRR(3)*100;
        
        C1 = SpectralClustering(Q1,numClust);
        result_LatLRR1 = ClusteringMeasure(truthF, C1);
        %         C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        %         result_LatLRR = EvaluationMetrics(truthF, C);  % res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];
        %c=1;
        acc1(iter_c)    = result_LatLRR1(1)*100;
        nmi1(iter_c) = result_LatLRR1(2)*100;
        purity1(iter_c)= result_LatLRR1(3)*100;
        
        C2 = SpectralClustering(Q2,numClust);
        result_LatLRR2 = ClusteringMeasure(truthF, C2);
        %         C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        %         result_LatLRR = EvaluationMetrics(truthF, C);  % res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];
        %c=1;
        acc2(iter_c)    = result_LatLRR2(1)*100;
        nmi2(iter_c) = result_LatLRR2(2)*100;
        purity2(iter_c)= result_LatLRR2(3)*100;
        
    end
    mean_ACC = mean(acc);           results(1) = mean_ACC;
    mean_NMI = mean(nmi);           results(2) = mean_NMI;
    mean_Purity = mean(purity);     results(3) = mean_Purity; 
    
    mean_ACC1 = mean(acc);           results1(1) = mean_ACC1;
    mean_NMI1 = mean(nmi);           results1(2) = mean_NMI1;
    mean_Purity1 = mean(purity);     results1(3) = mean_Purity1; 
    
    mean_ACC1 = mean(acc);           results2(1) = mean_ACC1;
    mean_NMI1 = mean(nmi);           results2(2) = mean_NMI1;
    mean_Purity1 = mean(purity);     results2(3) = mean_Purity1; 
    
    fprintf(' %5.4f \t  %5.4f \t   %5.4f \t   \n',mean_ACC, mean_NMI,mean_Purity); 
    fprintf(' %5.4f \t  %5.4f \t   %5.4f \t   \n',mean_ACC1, mean_NMI1,mean_Purity1); 
    fprintf(' %5.4f \t  %5.4f \t   %5.4f \t   \n',mean_ACC1, mean_NMI1,mean_Purity1); 
    
        
    input=[input,results(1),results(2),results(3)];

else
    for iter_c=1:T %
        %         C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        %         result_CLU = ClusteringMeasure(truthF, C);
        C = SpectralClustering(Q,numClust);% C = kmeans(U,numClust,'EmptyAction','drop');  % (165*1)
        result_LatLRR = EvaluationMetrics(truthF, C);  % res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];
        %c=1;
        acc(iter_c)   = result_LatLRR(1)*100;
        nmi(iter_c)   = result_LatLRR(2)*100;
        purity(iter_c)= result_LatLRR(3)*100;
        Fscore(iter_c)= result_LatLRR(4)*100;
        Preci(iter_c) = result_LatLRR(5)*100;
        Recall(iter_c)= result_LatLRR(6)*100;
        ARI(iter_c)   = result_LatLRR(7)*100;
        
        C1 = SpectralClustering(Q1,numClust);                       C2 = SpectralClustering(Q2,numClust);
        result_LatLRR1 = EvaluationMetrics(truthF, C1);             result_LatLRR2 = EvaluationMetrics(truthF, C2);  
        %c=1;
        acc1(iter_c)   = result_LatLRR1(1)*100;                     acc2(iter_c)   = result_LatLRR2(1)*100;
        nmi1(iter_c)   = result_LatLRR1(2)*100;                     nmi2(iter_c)   = result_LatLRR2(2)*100;
        purity1(iter_c)= result_LatLRR1(3)*100;                     purity2(iter_c)= result_LatLRR2(3)*100;
        Fscore1(iter_c)= result_LatLRR1(4)*100;                     Fscore2(iter_c)= result_LatLRR2(4)*100;
        Preci1(iter_c) = result_LatLRR1(5)*100;                     Preci2(iter_c) = result_LatLRR2(5)*100;
        Recall1(iter_c)= result_LatLRR1(6)*100;                     Recall2(iter_c)= result_LatLRR2(6)*100;
        ARI1(iter_c)   = result_LatLRR1(7)*100;                     ARI2(iter_c)   = result_LatLRR2(7)*100;
    end
    mean_ACC    = mean(acc);          results(1) = mean_ACC;  
    mean_NMI    = mean(nmi);          results(2) = mean_NMI;  
    mean_Purity = mean(purity);       results(3) = mean_Purity; 
    mean_Fscore = mean(Fscore);       results(4) = mean_Fscore;  
    mean_Preci  = mean(Preci);        results(5) = mean_Preci; 
    mean_Recall = mean(Recall);       results(6) = mean_Recall;  
    mean_ARI    = mean(ARI);          results(7) = mean_ARI;    
    
    mean_ACC1    = mean(acc1);          results1(1) = mean_ACC1;        mean_ACC2    = mean(acc2);          results2(1) = mean_ACC2;
    mean_NMI1    = mean(nmi1);          results1(2) = mean_NMI1;        mean_NMI2    = mean(nmi2);          results2(2) = mean_NMI2;
    mean_Purity1 = mean(purity1);       results1(3) = mean_Purity1;     mean_Purity2 = mean(purity2);       results2(3) = mean_Purity2; 
    mean_Fscore1 = mean(Fscore1);       results1(4) = mean_Fscore1;     mean_Fscore2 = mean(Fscore2);       results2(4) = mean_Fscore2; 
    mean_Preci1  = mean(Preci1);        results1(5) = mean_Preci1;      mean_Preci2  = mean(Preci2);        results2(5) = mean_Preci2; 
    mean_Recall1 = mean(Recall1);       results1(6) = mean_Recall1;     mean_Recall2 = mean(Recall2);       results2(6) = mean_Recall2; 
    mean_ARI1    = mean(ARI1);          results1(7) = mean_ARI1;        mean_ARI2    = mean(ARI2);          results2(7) = mean_ARI2; 
    fprintf(' Lv+Lc: %5.4f \t  %5.4f \t   %5.4f \t    %5.4f \t  %5.4f \t   %5.4f \t  %5.4f \t  \n',mean_ACC, mean_NMI,mean_Purity,mean_Fscore,mean_Preci,mean_Recall,mean_ARI);
    fprintf(' Lv:    %5.4f \t  %5.4f \t   %5.4f \t    %5.4f \t  %5.4f \t   %5.4f \t  %5.4f \t  \n',mean_ACC1, mean_NMI1,mean_Purity1,mean_Fscore1,mean_Preci1,mean_Recall1,mean_ARI1);
    fprintf(' Lc:    %5.4f \t  %5.4f \t   %5.4f \t    %5.4f \t  %5.4f \t   %5.4f \t  %5.4f \t  \n',mean_ACC2, mean_NMI2,mean_Purity2,mean_Fscore2,mean_Preci2,mean_Recall2,mean_ARI2);
    
    input=[input,results(1),results(2),results(3),results(4),results(5),results(6),results(7)];
end