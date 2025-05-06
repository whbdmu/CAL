function [a] = full_tr(tr,flag)
%              full_tr(tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX),  beta*ones(length(rX),1)))
% 输入: 
% tr=tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX), beta*ones(length(rX),1))=t
% flag

% Y=reshape(full_tr(tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX),beta*ones(length(rX),1))),[N,N,V]);  

%Converts TR-tensor to multi array
%
%
%
%---------------------------
% node:G(i),  d=5, n=rx=[11,15,11,15,3], r=[8,8,8,8,8]=bet=R
node=tr.node; d=tr.d; n=tr.n;   r=tr.r;



if 1        
    [~,id]=min(r); % 最小r值得索引 id=1
    
    
    a=node{id};
    for i=[id+1:d,1:id-1]
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    a=reshape(a,[r(id),prod(n),r(id)]);
    a=permute(a,[2,3,1]);
    a=reshape(a,[prod(n),r(id)*r(id)]);
    temp=eye(r(id),r(id));
    a=a*temp(:);
    
    a=reshape(a,n([id:d,1:id-1])');
    a=shiftdim(a,d-id+1);
    a=a(:);
    
    
end


if 0    
    a=node{1};
    for i=2:d
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    a=reshape(a,[r(1),prod(n),r(1)]);
    a=permute(a,[2,3,1]);
    a=reshape(a,[prod(n),r(1)*r(1)]);
    temp=eye(r(1),r(1));
    a=a*temp(:);
    
end

if 0
    a=node{2};
    for i=3:d
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    
    a=reshape(a,[r(2),prod(n(2:end)),r(1)]);
    a=permute(a,[1,3,2]);
    a=reshape(a,[r(2)*r(1), prod(n(2:end))]);
    
    b=node{1};
    b=permute(b,[2,3,1]);
    b=reshape(b,[n(1),r(2)*r(1)]);
    
    a=b*a;
    a=a(:);
    
    
    
end






if (nargin>1)&& (flag==1)
    a = reshape(a, n');
end;



return;

% 根据给定的代码片段，这段代码是一个矩阵张量网络（Matrix Product State, MPS）的运算过程。
% 根据不同的条件，它执行了不同的计算步骤。
% 
% 首先，根据 `if 1` 的条件，代码执行了一系列矩阵乘法和形状变换操作。具体步骤如下：
% 
% 1. 找到具有最小 `r` 值的索引 `id`。
% 2. 获取 `node{id}` 对应的矩阵 `a`。
% 3. 对于 `i` 从 `id+1` 到 `d`，以及从 1 到 `id-1` 的循环：
%    - 重塑矩阵 `cr` 的形状，以便进行矩阵乘法。如果 `i` 等于 `d`，则重塑为 `[r(i),n(i)*r(1)]` 形状，
%      否则重塑为 `[r(i),n(i)*r(i+1)]` 形状。
%    - 将矩阵 `a` 重塑为 `[numel(a)/r(i),r(i)]` 形状。
%    - 将矩阵 `a` 与矩阵 `cr` 相乘，更新 `a`。
% 4. 将 `a` 重塑为 `[r(id),prod(n),r(id)]` 形状。
% 5. 使用 `permute` 函数将维度重新排列为 `[2,3,1]`。
% 6. 将 `a` 重塑为 `[prod(n),r(id)*r(id)]` 形状。
% 7. 创建一个 `r(id)` × `r(id)` 的单位矩阵 `temp`。
% 8. 将 `temp` 展开为列向量，并将其与 `a` 相乘，更新 `a`。
% 9. 将 `a` 重塑为 `n([id:d,1:id-1])'` 形状。
% 10. 使用 `shiftdim` 函数将维度移动到正确的位置。
% 11. 最后，将 `a` 展平为列向量。
% 
% 接下来，根据 `if 0` 的条件，代码执行了另外两个计算步骤。具体步骤如下：
% 
% 1. 第一个 `if 0` 部分，类似于前面的步骤，但是从 `node{1}` 开始，然后依次处理每个节点。
%    最后，将结果 `a` 重塑为 `[r(1),prod(n),r(1)]` 形状，并展平为列向量。
% 2. 第二个 `if 0` 部分，类似于前面的步骤，但是从 `node{2}` 开始，然后依次处理每个节点。
%    最后，将结果 `a` 重塑为 `[r(2)*r(1), prod(n(2:end))]` 形状，并与 `node{1}` 相乘。
% 
% 最后，根据 `if (nargin>1)&& (flag==1)` 的条件，如果输入参数中包含 `flag` 并且其值为 1，
% 代码将根据指定的维度 `n` 对结果 `a` 进行形状重塑。
% 
% 整个代码段的目的是对一个张量网络进行计算，并返回结果 `a`。
% 根据具体的输入数据和条件，最终结果的维度可能会有所不同。
