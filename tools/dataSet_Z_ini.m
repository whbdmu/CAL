function [X,Z_ini] = dataSet_Z_ini(X,k,ind_folds)
for iv = 1:length(X) % 视图数
    X1 = X{iv}; % (维度*样本数)
    X1 = NormalizeFea(X1,0);               % 归一化，按列——样本
    ind_0 = find(ind_folds(:,iv) == 0);    % (N*MR(m)*1)  % 第iv个视图中 数值为0的索引——缺失实例的索引
    ind_1 = find(ind_folds(:,iv) == 1);    % (1032*1)     % 第iv个视图中 数值为1的索引——可用实例的索引
    % length(ind_0)+length(ind_1)=samples=1474
    %% 缺失视角补0
    X1(:,ind_0) = 0; % (48*1474) X1中索引为ind_0中的那一列都为0，列对应的是样本
    Y{iv} = X1;      % 实例集——缺失实例所在列为0    % 一列一个样本
    
    %% ---------- 初始KNN图构建------------ %
    X1(:,ind_0) = []; % X1(48*1474),删除X1中的缺失实例，得 可用实例集(48*1032)
    options = [];
    options.NeighborMode = 'KNN';
    options.k = k; % 10
    options.WeightMode = 'Binary';      % Binary  HeatKernel (0,1)
    Z1 = full(constructW(X1',options)); % 稀疏矩阵-->全矩阵  (可用*可用)  % X1' (1032*48)
    linshi_G = diag(ind_folds(:,iv)); % (1474*1474)(n*n)
    linshi_G(:,ind_0) = []; %  缺失实例所在列删除 (1474*1032)(n*可用)
    Z_ini{iv} = linshi_G*max(Z1,Z1')*linshi_G'; % (n*可用)(可用*可用)(可用*n)
    
    clear Z1 linshi_G
end  % for iv = 1:length(X)
%clear X X1 ind_0
X = Y;  % X中缺失实例所在列 元素都为0 (样本数*维度)
clear Y
end


