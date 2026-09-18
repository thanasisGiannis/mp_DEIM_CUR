function [U,S,V] = rsvd_half(A,k,p,q)
%RSVD_DOUBLE  Rank-k randomized SVD in double precision.
%
% [U,S,V] = rsvd_half(A,k)
% [U,S,V] = rsvd_half(A,k,p,q)
%
% A : m-by-n matrix
% k : requested rank
% p : oversampling parameter (default 10)
% q : number of power iterations (default 0)

options.format = 'h';

if nargin < 3 || isempty(p), p = 10; end
if nargin < 4 || isempty(q), q = 0;  end

A = chop(A,options);
[m,n] = size(A);
ell = min(k + p, min(m,n));

% Random sampling.
Omega = chop(randn(n, ell),options);
Y = chop(A * Omega,options);

% Optional subspace iteration, with reorthogonalization.
for j = 1:q
    [Q,~] = qr(Y, 0); Q = chop(Q,options);
    Y = chop(A' * Q,options);
    [Q,~] = qr(Y, 0); Q = chop(Q,options);
    Y = chop(A * Q,options);
end

% Orthonormal basis for the sampled range.
[Q,~] = qr(Y, 0); Q = chop(Q,options);

% Small SVD and truncation.
B = chop(Q' * A,options);
[Uh,S,V] = svd(B, 'econ');
Uh = chop(Uh,options);
S  = chop(S,options);
V  = chop(V,options);

U = chop(Q * Uh(:,1:k),options);
S = S(1:k,1:k);
V = V(:,1:k);
end
