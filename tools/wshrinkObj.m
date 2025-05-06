function [x,objV] = wshrinkObj(x,           rho,         sX,   isWeight,mode)
%    [Pv, objV] = wshrinkObj(Zv + 1/miu*Av,1/miu,[Nsamp,Nsamp,nv],0,      1);
% Zv(13036056*1),Av(13036056*1)

if isWeight == 1
    % ��� isWeight ���� 1�������Ȩ�� C��C ��ֵ�Ǹ���������� sX ��ά�ȼ���õ���
%     C = 2*sqrt(2)*sqrt(sX(3)*sX(2));
    C = sqrt(sX(3)*sX(2)); % (nv*Nsamp)^(1/2)
end
% ���û��ָ�� mode ������Ĭ����Ϊ 1��mode ��������ѡ��ͬ����Ƭ����
if ~exist('mode','var')
    % mode = 1�ǲ���lateral slice�ķ���
    % mode = 2�ǲ���front slice�ķ���
    % mode = 3�ǲ���top slice�ķ���
    mode = 1;
end

% ���������� x ����ת��Ϊ��СΪ sX ������ X
X=reshape(x,sX); % (1474*1474*6)
%% ִ����Ƭ
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);
else
    Y = X;
end
%% 

% �� Y ���п��ٸ���Ҷ�任���õ� Yhat
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
    % ���� endValue ��ֵ��ȡ n3/2+1 ���������֣���ѭ���������� Yhat �е�ÿ����Ƭ
    endValue = int16(n3/2+1);     % ����2 Ϊʲô�������أ�
    for i = 1:endValue
        % �Ե�ǰ��Ƭ��������ֵ�ֽ⣨SVD�����õ��������������� uhat������ֵ���� shat ���������������� vhat
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        
        if isWeight
            % ��� isWeight Ϊ�棬�����Ȩ�� weight�������� tau �� shat ��������ֵ����
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat = max(shat - tau,0);
        end                 
        
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat'; % Yhat=Z+A/miu��Eq(13)��
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
%             conj�󸴹����
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
 