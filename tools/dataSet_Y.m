function [X,X_E,G,Y,Omega,omega] = dataSet_Y(X,V,ind_folds)

omega = cell(1,V);
for iv = 1:length(X) % 视图数
    X1 = X{iv}; % (维度*样本数)
    X1 = NormalizeFea(X1,0);               % 归一化，按列——样本
    ind_0 = find(ind_folds(:,iv) == 0);    % (N*MR(m)*1)  % 第iv个视图中 数值为0的索引——缺失实例的索引
    ind_1 = find(ind_folds(:,iv) == 1);    % (1032*1)     % 第iv个视图中 数值为1的索引——可用实例的索引
    % 缺失视角补0
    X1(:,ind_0) = 0; % (48*1474) X1中索引为ind_0中的那一列都为0，列对应的是样本
    X1_0{iv} = X1;      % 实例集——缺失实例所在列为0    % 一列一个样本
    X2 = X1;
    % 删除缺失实例
    X2(:,ind_0) = [];
    X_E{iv} = X2;  % 可用实例集 (dv*可用)
    % 观测子集矩阵
    omega{iv} = zeros(size(X1));
    view_v_complete_instances = setdiff(1:size(X{iv},2), ind_0);
    omega{iv}(:,view_v_complete_instances) = 1;
    %         X_E2{iv} = omega{iv}.*X{iv};
    %         X_E2{iv} = NormalizeFea(X_E2{iv},0);
    %         if isequal(X1_0{iv},X_E2{iv})   fprintf('相等');
    %         else fprintf('不相等');   end
    
    linshi_G = diag(ind_folds(:,iv)); % (1474*1474)(n*n)
    linshi_G(:,ind_0) = [];
    G{iv} = linshi_G;
end  % for iv = 1:length(X)
%clear X X1 ind_0
X = X1_0;  % X中缺失实例所在列 元素都为0 (样本数*维度)
Y = X{1};
Omega = omega{1};
for v = 2:V
    Y = [Y;X{v}];
    Omega = [Omega;omega{v}];
end

end