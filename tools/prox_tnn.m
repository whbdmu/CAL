function [X, tnn, trank] = prox_tnn(M,        rho,  p,mode)
%        [M,  ~,    ~]   = prox_tnn(Z+Q1/mu,beta/mu,p,mode);
% beta=ones(n1,1)， p=0.6， mode=2， r=1.3  mu=1e-4;
% M(400*400*4), rho(400*1)10000, p=0.6, modde=2  

%this function is used to update E of our model,E is the tensor

% The proximal operator of the tensor nuclear norm of a 3 way tensor
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2  Eq(13)
%
% Y     -    n1*n2*n3 tensor
%
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X
%

% 
% 切片方式
if mode == 1
    Y=X2Yi(M,3);
elseif mode == 3
    Y=shiftdim(M, 1);
else  % mode=2
    Y = M; % M=Z+Q1/mu (400*400*4)
end

[n1,n2,n3] = size(Y); % 400 400 4
n12 = min(n1,n2); % 400

Y = fft(Y,[],3); % (400*400*4)

U = zeros(n1,n12,n3); % (400*400*4)
V = zeros(n2,n12,n3); % (400*400*4)
S = zeros(n12,n12,n3); % (400*400*4)
trank = 0;
for i = 1 : n3  % 1:4
    [U(:,:,i),s,V(:,:,i)] = svd(Y(:,:,i),'econ'); % econ 经济型奇异值分解
    % U(400*0*4) s(400*1) V(400*0*4)
    % i=1,s=0
    s = diag(s);
    % rho(400*1)10000 p=0.6
    s = solve_Lp_w(s, rho, p); % 400*1  
    S(:,:,i) = diag(s); % (0*0*4)
    tranki = length(find(s~=0)); % s不等于的个数
    trank = max(tranki,trank);
end
U = U(:,1:trank,:); % (400*0*4)
V = V(:,1:trank,:);
S = S(1:trank,1:trank,:);

U = ifft(U,[],3);
S = ifft(S,[],3);
V = ifft(V,[],3);

X = tprod( tprod(U,S), tran(V)); % 张量的积  张量的转置

S = S(:,:,1);
tnn = sum(S(:)); % return the tensor nuclear norm of X
