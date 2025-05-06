function [X1, O1, X2, O2] = DataPreparing(data, ind_folds)
% data(1*3)(样本*维度)(3560/3631/3068*416)
% index (1*3)(352/302/294*1)

%        X{k} = X1{k}*O1{k} + X2{k}*O2{k}  
% Eq(14) X(k) = X(k)o*P(k)o + X(k)u*P(k)u
% X(k)o == X1{k} 可用样本(D,N_1)，
% X(k)u == X2{k} 缺失样本(D,N_i)

% [D, N] = size(X{k}), [D, N_1] = size(X1{k}), [D, N_i] = size(X2{k}),
% [N_1, N] = size(O1{k}), [N_i, N] = size(O2{k})
% N = N_1 + N_i

K = length(data); %numOfView           3
N = size(data{1}, 2); %numOfSample    416  data{1} 3562*416
X1 = cell(K,1); %the complete parts   (3*1)   可用
X2 = cell(K,1); %the missing parts    (3*1)   缺失
% 两个 置换矩阵
O1 = cell(K, 1); %                    (3*1)   
O2 = cell(K, 1); %                    (3*1)
Xc = cell(K,1);  %                    (3*1)

index=cell(1,K);
for i = 1:K
   index{i} = find(ind_folds(:,i) == 1); 
end

for k = 1:K  % 3
    data{k} = data{k};  % 视图1 2 3 (3560/3631/3068*416)

    W1 = ones(N, 1);  % (样本416*1) 全1
    W1(index{k}, 1) = 0;  % 缺失实例 索引为1；可用实例 索引为0
    
    ind_1 = W1 == 1; % (416*1) 缺失实例-1（7、16），可用-0（1 2 3 4 5 6）
    W2 = eye(N);  % (416*416)的单位矩阵
    W2(ind_1, :) = []; % (352*416)(N_1*N) 可用实例数N_1 = 352
    O1{k} = W2; % 可用实例的投影矩阵(N_1*N)  Po{k}
    
    W3 = zeros(N, 1); % (样本416*1) 全0
    W3(index{k}, 1) = 1; % 可用实例 索引为1；缺失实例 索引为0
    
    
    ind_2 = W3 == 1; % (416*1) 缺失实例-1（1 2 3 4 5 6），可用-0（7 16）
    W4 = eye(N);
    W4(ind_2, :) = []; % (64*416)(N_i*N)
    O2{k} = W4; % 缺失实例的投影矩阵
    
    data{k} = double(data{k}); % 双精度
    data{k}(isnan(data{k})) = 0; % nan值 被置为0
    Xc{k} = data{k} * O1{k}'; % (3*1)(维度*可用样本)
    [X1{k}] = NormalizeData(Xc{k});  % 归一化
    fillV = repmat(mean(X1{k}, 2), 1, size(O2{k}, 1));
    % O2{k} (缺失样本*416)(N_i*N)   fillV (维度*缺失样本)
    % 计算 X1{k} 按列求平均值，并使用 repmat 函数将平均值扩展为与 O2{k} 大小相同的矩阵，并将结果存储在 fillV 中。
    [X2{k}] = NormalizeData(fillV);
end
end
% X1 (3*1)(维度*可用样本)(D*N_1)
% X2 (3*1)(维度*缺失样本)(D*N_i)