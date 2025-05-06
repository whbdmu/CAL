
function [t] = tensor_ring_als(c,Rank,mu)
%              tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX),beta*ones(length(rX),1))
% 输入: c=reshape(Z_tensor-1/rho3*L_tensor,rX), 张量Z/L_tensor(N*N*V)-->(5*1)5D
%       Rank=beta*ones(length(rX),1)
%       mu=

% Y=reshape(full_tr(tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX),beta*ones(length(rX),1))),[N,N,V]);   

% Tensor Ring toolbox
% written by Qibin Zhao, BSI, RIKEN

if (nargin == 0)
    t.d    = 0;
    t.r    = 0;
    t.n    = 0;
    t.node = 0;                    % empty tensor
    t = class(t, 'tensor_ring');
  
    return;
end
% ip = inputParser;
% ip.addParamValue('Tol', 1e-6, @isscalar);
% ip.addParamValue('Rank', [], @ismatrix);
% ip.addParamValue('MaxIter', 1000, @isscalar);
% ip.addParamValue('weight', 0.01, @isscalar);
% ip.addParamValue('lambda', 0.01, @isscalar);
%  ip.parse(varargin{2:end});
% 
% Tol = ip.Results.Tol;
% Rank = ip.Results.Rank;
% MaxIter = ip.Results.MaxIter;
% lambda=ip.Results.lambda;


 maxit=8;
 Tol=1e-5;
 
%  Rank=[8,8,8,8,10];
        n = size(c); % c:5-D double 五维数组，n=rx=[11,15,11,15,3]； val(:,:,1,1,1)~val(:,:,11,15,3)
        n = n(:); % 转为列向量
        d = numel(n); % 5 计算元素总和
        node=cell(1,d);
        r=Rank(:); % ；列向量，元素均相等 为 beta=j=R=8;
        %% TR分解： r(i):TR秩， n=rx=I，  G(i)=Z(i)=(Ri*Ii*Ri+1)
        for i=1:d-1                          
            % 在每次循环中，使用 randn 函数生成一个随机矩阵，然后将其赋值给 node 数组的第 i 个元素。
            node{i}=randn(r(i),n(i),r(i+1)); % 生成的随机矩阵的大小由 r(i)、n(i) 和 r(i+1) 决定。  (8*  *8)
        end                             % R1 I1 R2, R2 I2 R3, R3 I3 R4, R4 I4 R5, 
        node{d}=randn(r(d),n(d),r(1));  % R5 I5 R1
        od=[1:d]'; % 列向量 从1-5
        err=1;
        for it=1:maxit % 1:8
            err0=err; % 1
            if it>1
                c=shiftdim(c,1); % 维度左移1位，val(:,:,15,3,11)
                od=circshift(od,-1); % 循环向左移位1个位置，列向量（2,3,4,5,1）
            end % it=1
            n_od_1 = n(od(1)); % 11
            num_c = numel(c); % 81675
            nn = numel(c)/n(od(1)); % 7425=81675/11
            c=reshape(c,n(od(1)),numel(c)/n(od(1)));% (11*7425) (15*7425)
            b=node{od(2)};  % 2: (8*15*8)
            for k=3:d
                j=od(k); % j=k; k+1
                br=node{j}; % 3: (8* 11/15/3 *8); (8* 15/3/11 *8)
                br=reshape(br,[r(j),numel(br)/r(j)]); % [8,88=8*11*8/8]  % numel(br)  [8,120=8*15]   [8,24=8*3]
                b=reshape(b,[numel(b)/r(j),r(j)]);    % [120=8*15*8/8,8] % numel(b)   [1320=120*88/8,8]  [19800=1320*120/8,*8]
                b=b*br; % (120*88)  (1320*120) (19800*24)
            end % prod()所有元素计算得到的乘积
            b=reshape(b,[r(od(2)),prod(n(od(2:end))),r(od(1))]); % [8,15*11*15*3=7425,8]
            b=permute(b,[1,3,2]); % 矩阵 b 进行维度重排,[1,3,2]指定的目标维度顺序，(8*8*7425)
            b=reshape(b,[r(od(2))*r(od(1)), prod(n(od(2:end)))]); % [8*8=64,7425]
%              a=c*b'*inv(b*b'+mu*eye(r(od(2))*r(od(1))));
               a=c/b; % (11*64)=(11*7425)/(64*7425)
              % Un(i,:)=Un(i,:).*(Xn(i,idx)*Uneq(idx,:))./(Un(i,:)*(Uneq(idx,:)'*Uneq(idx,:))); 

%             inv(Uneq(idx,:)'*Uneq(idx,:)+para.mu*eye(R(n)*R(n+1)))
            
            err=norm(c-a*b,'fro')/norm(c(:));
            a=reshape(a,[n(od(1)),r(od(2)),r(od(1))]); % (11*8*8)
            node{od(1)}=permute(a,[3,1,2]);
            s=norm(node{od(1)}(:));
            node{od(1)}=node{od(1)}./s;
            
%             fprintf('it:%d, err=%f\n',it,err);
%             error(it)=err;
            if abs(err0-err)<=1e-5 && it>=2*d && err<=Tol
                break;
            end
            c=reshape(c,n(od)');  % 5-D double
            
        end

        node{od(1)}=node{od(1)}.*s;
        t.node=node;
        t.d=d;
        t.n=n;
        t.r=r;    
        return;
    end