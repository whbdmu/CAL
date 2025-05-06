function [x,objV] = wshrinkObj(x,           rho,         sX,   isWeight,mode)
%    [Pv, objV] = wshrinkObj(Zv + 1/miu*Av,1/miu,[Nsamp,Nsamp,nv],0,      1);
% Zv(13036056*1),Av(13036056*1)

if isWeight == 1
    % 如果 isWeight 等于 1，则计算权重 C。C 的值是根据输入参数 sX 的维度计算得到的
%     C = 2*sqrt(2)*sqrt(sX(3)*sX(2));
    C = sqrt(sX(3)*sX(2)); % (nv*Nsamp)^(1/2)
end
% 如果没有指定 mode 参数，默认设为 1。mode 参数用于选择不同的切片方法
if ~exist('mode','var')
    % mode = 1是采用lateral slice的方法
    % mode = 2是采用front slice的方法
    % mode = 3是采用top slice的方法
    mode = 1;
end

% 将输入向量 x 重新转换为大小为 sX 的数组 X
X=reshape(x,sX); % (1474*1474*6)
%% 执行切片
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);
else
    Y = X;
end
%% 

% 对 Y 进行快速傅里叶变换，得到 Yhat
Yhat = fft(Y,[],3);
% weight = C./sTrueV+eps;
% weight = 1;
% tau = rho*weight;
objV = 0;
if mode == 1     % lateral slice
    n3 = sX(2);
elseif mode == 3 % top slice
    n3 = sX(1);
else             % front slice
    n3 = sX(3);
end

if isinteger(n3/2)
    % 根据 endValue 的值（取 n3/2+1 的整数部分），循环迭代处理 Yhat 中的每个切片
    endValue = int16(n3/2+1);     % 问题2 为什么不跑完呢？
    for i = 1:endValue
        % 对当前切片进行奇异值分解（SVD），得到左奇异向量矩阵 uhat、奇异值矩阵 shat 和右奇异向量矩阵 vhat
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        
        if isWeight
            % 如果 isWeight 为真，则计算权重 weight，并根据 tau 对 shat 进行软阈值处理
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat = max(shat - tau,0);
        end                 
        
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat'; % Yhat=Z+A/miu，Eq(13)下
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
            objV = objV + sum(shat(:));
        end
    end
    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');
    if isWeight
       weight = C./(diag(shat) + eps);
       tau = rho*weight;
       shat = soft(shat,diag(tau));
    else
       tau = rho;
       shat = max(shat - tau,0);
    end
    
    objV = objV + sum(shat(:));
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
else
   endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        if isWeight
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat = max(shat - tau,0);
        end
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
%             conj求复共轭函数
            objV = objV + sum(shat(:));
        end
    end 
end

Y = ifft(Yhat,[],3);
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X = shiftdim(Y, 2);
else
    X = Y;
end

x = X(:);

end
 