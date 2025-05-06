function X = prox_l21(B,lambda)
%     E{i} = prox_l21(X{i}-X{i}*Z{i} + F1{i}/rho, 1/rho);

% The proximal operator of the l21 norm of a matrix
% l21 norm is the sum of the l2 norm of all columns of a matrix 
%
% min_X lambda*||X||_{2,1}+0.5*||X-B||_2^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
%

X = zeros(size(B));
for i = 1 : size(X,2)
    nxi = norm(B(:,i));
    if nxi > lambda  
        X(:,i) = (1-lambda/nxi)*B(:,i);
    end
end

% 这段代码实现了矩阵的l21范数的近似最小化（proximal operator），具体解决的问题是：
% 
% min_X lambda*||X||_{2,1}+0.5*||X-B||_2^2
% 
% 下面是代码的详细解释：
% 
% 首先，初始化解矩阵 X 为与输入矩阵 B 相同尺寸的零矩阵。
% 
% 然后，通过循环遍历 X 的每一列（特征向量），计算每列的l2范数。使用 norm() 函数计算矩阵 B 的每一列的l2范数，得到 nxi。
% 
% 对于每一列，判断其l2范数是否大于 lambda。如果是，则进行近似最小化操作。
% 
% 近似最小化操作的具体步骤是，计算新的列向量 X(:,i)。如果 l2范数 nxi 大于 lambda，
% 则根据公式 (1-lambda/nxi)*B(:,i) 计算新的列向量，并将其赋值给 X(:,i)。这个公式通过缩放原始列向量，使其满足l2范数约束。
% 
% 最后，返回求解得到的矩阵 X，表示满足约束条件的近似最小化解。
% 
% 通过执行这段代码，可以求解矩阵的l21范数的近似最小化问题，得到满足约束条件的解矩阵 X。