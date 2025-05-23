function x = prox_l1(b,lambda)
% S{i} = prox_l1(O{i}*Z{i}*O{i}'-B{i}+F2{i}/rho, theta/rho); % F2=F1

% The proximal operator of the l1 norm
% 
% min_x lambda*||x||_1+0.5*||x-b||_2^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

x = max(0,b-lambda)+min(0,b+lambda);