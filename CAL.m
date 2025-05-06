function [Z,L,Lc,obj,iter] = CAL(Xo,Po,Xu,Pu,numClust,truthF,paras)

V = length(Xo); 
N = size(Po{1},2); % Sample number
c = numClust;
% Initialization
lambda1 = paras.lambda1;
lambda2 = paras.lambda2;
p = paras.p;
mode = paras.mode;
max_iter = paras.max_iter;

rho = paras.rho;
miu = paras.miu;
yita = paras.yita;
tol = paras.tol;


Z = cell(1,V);
X = cell(1,V);
Xc = cell(1,V);
for v = 1:length(X)
    Z{v} = zeros(N, N);
    X{v} = Xo{v} * Po{v};  % (dim*N)
    Xc{v} = Xo{v} * Po{v};
    B1{v} = zeros(size(X{v})); % X{v}-X{v}Z{v}-B1{v};   
    L{v} = zeros(size(Z{v}));  % Z{v}=L{v}+S{v};  
    S{v} = zeros(size(Z{v}));  
    B2{v} = zeros(size(L{v})); % L{v}=Lc+B{v};     
    A1{v} = zeros(size(L{v}));
    % Lagrange multipliers
    Q1{v} = zeros(size(Z{v})); % Z{v}-L{v}-S{v};
    Q2{v} = zeros(size(L{v})); % L{v}-Lc-B{v};   
    Y1{v} = zeros(size(L{v}));
end
Lc = zeros(size(L{v}));
Ac = Lc;  
Yc = Lc;
dims=[N,N,(V+1)];
LL_tensor = cat(3, L{:,:},Lc); % (N,N,(V+1))
beta1 = ones(1,(V+1))';
sX = [N,N,(V+1)];

% Optimization 

for iter = 1:max_iter
    X_pre = X;
    Z_pre = Z;
    L_pre = L;
    Lc_pre = Lc;

    % -------- Xu(v)---------%                                                    
    M = cell(1,V);
    for v = 1 : V
        M{v} = (Z{v} - eye(N)) * (Z{v} - eye(N))';
        tmp1 = Pu{v} * M{v} * Pu{v}';
        tmp2 = -(Xo{v} * Po{v} * M{v} * Pu{v}' + B1{v}*(Z{v} - eye(N))'*Pu{v}');
        Xu{v} = tmp2 / tmp1;
        [Xu{v}] = NormalizeData(Xu{v});
        X{v} = Xc{v} + Xu{v} * Pu{v};
    end
    
    % --------- Z{v} ------%
    for v = 1:V
        tmpz1 = X{v}'*X{v}+rho*eye(size(X{v},2));
        tmpz2 = X{v}'*X{v}-X{v}'*B1{v}+rho*(L{v}+S{v})-Q1{v};
        linshi3=max(inv(tmpz1)*tmpz2,0);
        Z1 = zeros(size(linshi3));
        for is = 1:size(linshi3,1)
            ind_c = 1:size(linshi3,1);
            ind_c(is) = [];
            Z1(is,ind_c) = EProjSimplex_new(linshi3(is,ind_c));
        end
        Z{v} = Z1;
    end
    clear tmpz1 tmpz2 linshi3
    
    % -------- B{v}-- ------- %
    for v=1:V
        tmpB1 = X{v}-X{v}*Z{v};          tmp1 = lambda1;
        tmpB2 = L{v}-Lc+1/rho*Q2{v};     tmp2 = lambda1/rho;
        B1{v} = max(0,tmpB1-tmp1) + min(0,tmpB1+tmp1);
        B2{v} = max(0,tmpB2-tmp2) + min(0,tmpB2+tmp2);
    end
    clear tmpB1 tmpB2 tmp1 tmp2 

    % -------- L{v}--------%
    for v = 1:V        
        tmpl2 = 1/(2*rho+miu);
        tmpl1 = rho*(Z{v}-S{v}) + Q1{v} + rho*(Lc+B2{v}) - Q2{v} + miu*A1{v} - Y1{v};
        L{v} = tmpl2 * tmpl1;
    end    
    clear tmpl1 tmpl2

    % -------- S{v} -------%
    S_tensor = cat(3, S{:,:});
    Z_tensor = cat(3, Z{:,:});
    L_tensor = cat(3, L{:,:});
    Q1_tensor = cat(3, Q1{:,:}); % Update S_tensor
    for v = 1:V
        S_tensor(:,:,v) = prox_l21(Z_tensor(:,:,v)-L_tensor(:,:,v)+1/rho*Q1_tensor(:,:,v),lambda2/rho); % 
        S{v} = S_tensor(:,:,v);
    end
    
    % -------- Tensor A1 -------%
    LL_tensor = cat(3, L{:,:},Lc); 
    A1_tensor = cat(3,A1{:,:},Ac);
    Y1_tensor = cat(3,Y1{:,:},Yc);
    tmpa_t = LL_tensor + 1/miu*Y1_tensor;  % (n*n*v)
    tmpa = tmpa_t(:); % (n*n*v,1)
    [myj, ~] = wshrinkObj_weight_lp(tmpa, beta1./miu,sX, 0,3,p); % beta1 = ones(1,V)';  sX = [N,N,V]; mode=3;
    A1_tensor = reshape(myj, sX);
    for v = 1:V
       A1{v} = A1_tensor(:,:,v);
    end
    Ac = A1_tensor(:,:,V+1);
    
    % -------- Lc -------% 
    tmp1 = 1/(rho*V+miu);
    tmp2=0;
    for v = 1:V
        tmp2 = tmp2 + rho*(L{v}-B2{v}) + Q2{v};
    end
    tmp3 = tmp2 + miu*Ac - Yc;
    Lc = tmp1*tmp3;
    clear tmp1 tmp2 tmp3
    
    % ----------Lagrange multipliers------------ Q
    LL_tensor = cat(3, L{:,:},Lc); 
    A1_tensor = cat(3,A1{:,:},Ac);
    Y1_tensor = cat(3,Y1{:,:},Yc);
    Y1_tensor = Y1_tensor + miu*(LL_tensor-A1_tensor);
    for v = 1:V
        Q1{v} = Q1{v} + rho*(Z{v}-L{v}-S{v});
        Q2{v} = Q2{v} + rho*(L{v}-Lc-B2{v});
        Y1{v} = Y1_tensor(:,:,v);
    end
    Yc=Y1_tensor(:,:,V+1);

    % ------- penalty parameter -----------%
    rho = min(rho*yita, 1e10);
    miu = min(miu*yita, 1e10);
    
    res=0;
    diff_X = 0;
    diff_Z = 0;
    diff_L = 0;
    diff_Lc = 0;
    Rec1=0; rec1=0; 
    Rec2=0; rec2=0;
    Rec_error(iter) = 0;
    
    % check convergence
    for v = 1:V
        Rec_error(iter) = Rec_error(iter) + norm(X{v} - X{v} * Z{v} - B1{v}, 'fro') ^ 2;
        % Matching error
        rec1 = Z{v}-L{v}-S{v};          
        Rec1 = max(Rec1,max(abs(rec1(:))));
        rec2 = L{v}-Lc-B2{v};            
        Rec2 = max(Rec2,max(abs(rec2(:))));
        % Iterative error
        diff_X = max(diff_X,max(abs(X{v}(:)-X_pre{v}(:)))); 
        diff_Z = max(diff_Z,max(abs(Z{v}(:)-Z_pre{v}(:))));
        diff_L = max(diff_L,max(abs(L{v}(:)-L_pre{v}(:))));
    end
    diff_Lc = max(diff_Lc,max(abs(Lc(:)-Lc_pre(:))));
    History(iter) = Rec_error(iter);
    if iter >=2
        res = abs((History(iter) - History(iter - 1)) / History(iter - 1));
        err = min([res,Rec1,Rec2,diff_X,diff_Z,diff_L,diff_Lc]);
        obj(iter) = err;
        if ( err < tol && iter>10)
            iter;
            break;
        end
    end  
end % for iter = 1:max_iter 

% Clustering
KK=0; 
for i=1:length(X)
    KK = KK + (abs(L{i})+(abs(L{i}))')/2;
end
Q = KK/length(X) + (abs(Lc)+(abs(Lc))')/2;
for iter_c=1:10 %
    C = SpectralClustering(Q,numClust);
    result_LatLRR = EvaluationMetrics(truthF, C);  % res = [acc, nmi, Pu, Fscore, Precision, Recall, ARI];
    acc(iter_c)    = result_LatLRR(1)*100;
    nmi(iter_c) = result_LatLRR(2)*100;
    purity(iter_c)= result_LatLRR(3)*100;
    Fscore(iter_c)= result_LatLRR(4)*100;
    Preci(iter_c) = result_LatLRR(5)*100;
    Recall(iter_c)= result_LatLRR(6)*100;
    ARI(iter_c)   = result_LatLRR(7)*100;
end
mean_ACC = mean(acc);
mean_NMI = mean(nmi);
mean_Purity = mean(purity);
mean_Fscore = mean(Fscore);
mean_Preci = mean(Preci);
mean_Recall = mean(Recall);
mean_ARI = mean(ARI);
fprintf(' %5.4s \t  %5.6s \t   %5.6s \t    %5.6s \t  %5.6s \t   %5.6s \t  %5.6s \t   %5.6s\n','ACC', 'NMI','Purity','Fscore','Pre','Recall,','ARI','iter');
fprintf(' %5.4f \t  %5.4f \t   %5.4f \t    %5.4f \t  %5.4f \t   %5.4f \t  %5.4f \t   %d    \n',mean_ACC, mean_NMI,mean_Purity,mean_Fscore,mean_Preci,mean_Recall,mean_ARI,iter);
end % function
