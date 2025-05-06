
clear memory
clear all
addpath('tools');
addpath('datasets');
rand('seed',5867);

%  Load data  
MR=[0.1,0.3,0.5,0.7];  % missing rate
m=1;

Dataname = "bbcsport4vbigRnSp";
load(Dataname);  
datafile = strcat(Dataname, '_percentDel_',num2str(MR(m)),'.mat');
load(datafile);
ind_idx=3;
ind_folds = folds{ind_idx};
truthF=double(truth); 

N = length(truthF);
V = length(X);
numClust = length(unique(truthF));
    
for lambda1 = 0.01
    for lambda2 = 0.01
        for p = 0.7
            paras.lambda1 = lambda1;
            paras.lambda2 = lambda2;
            paras.p = p;
            paras.rho = 0.01;
            paras.miu = 0.01;   
            paras.yita = 3;    
            paras.mode = 2;
            paras.max_iter = 100;
            paras.tol = 1e-6;
            fprintf('Datanmae:%s,   miss：%.4f,  λ1：%f,  λ2：%f,   rho：%f,   miu：%f,   yita：%f,    p：%f,    fold:%d   \n',Dataname,MR(m),paras.lambda1,paras.lambda2,paras.rho,paras.miu,paras.yita,paras.p,ind_idx);
            % Data preparation
            [Xo, Po, Xu, Pu] = DataPreparing(X, ind_folds);
            % Train
            t0=tic;
            [Z,L,Lc,obj,iter] = CAL(Xo,Po,Xu,Pu,numClust,truthF,paras);
            Time=toc(t0);
        end
    end 
end 





