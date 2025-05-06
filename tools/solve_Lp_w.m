function   x   =  solve_Lp_w( y, lambda, p )  % （奇异值，lambda/miu，p）
% Modified by Dr. xie yuan
% lambda here presents the weights vector
J     =   4;  %2
% tau is generalized thresholding vector 正则化阈值向量
% 正则化因子lambda和正则化指数p
tau   =  (2*lambda.*(1-p)).^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda).^((p-1)/(2-p));
x     =   zeros( size(y) ); % （v*1）
% i0 is the number of zeros after thresholding
i0    =   find( abs(y)>tau );

if length(i0)>=1
    % lambda  =   lambda(i0);
    y0    =   y(i0);
    t     =   abs(y0);
    lambda0 = lambda(i0);
    %lambda0 = lambda;  %lambda0 = lambda(i0);
    for  j  =  1 : J
        t    =  abs(y0) - p*lambda0.*(t).^(p-1);
    end
    x(i0)   =  sign(y0).*t;
end


% 函数实现：
% 
% 首先定义了一个常量J，它是阈值处理的迭代次数。然后计算了一个阈值向量tau，其中包括了正则化因子lambda和正则化指数p。
% 
% 接着，定义了一个全零向量x，其大小与输入向量y相同。
% 然后，通过寻找abs(y)大于阈值向量tau的元素的索引，计算了一个非零元素的索引向量i0。
% 
% 如果非零元素的索引向量i0的长度大于等于1，则将这些元素存储在向量y0中，并根据迭代次数J进行阈值处理，得到一个新的向量t。
% 最后，将处理后的结果存储在原向量中对应的位置i0内。
% 
% 总的来说，该函数实现了带权Lp正则化。
% 其主要思想是在输入向量中对非零元素进行阈值处理，使其满足正则化约束，并保持输入向量中的零元素不变。