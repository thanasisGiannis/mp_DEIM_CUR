function [U,S,V] = rsvd_double(A,k,p,q)
%RSVD_DOUBLE  Rank-k randomized SVD in double precision.
%
% [U,S,V] = rsvd_double(A,k)
% [U,S,V] = rsvd_double(A,k,p,q)
%
% A : m-by-n matrix
% k : requested rank
% p : oversampling parameter (default 10)
% q : number of power iterations (default 0)

if nargin < 3 || isempty(p), p = 10; end
if nargin < 4 || isempty(q), q = 0;  end

A = double(A);
[m,n] = size(A);
ell = min(k + p, min(m,n));

% Random sampling.
Omega = randn(n, ell);
Y = A * Omega;

% Optional subspace iteration, with reorthogonalization.
for j = 1:q
    [Q,~] = qr(Y, 0);
    Y = A' * Q;
    [Q,~] = qr(Y, 0);
    Y = A * Q;
end

% Orthonormal basis for the sampled range.
[Q,~] = qr(Y, 0);

% Small SVD and truncation.
B = Q' * A;
[Uh,S,V] = svd(B, 'econ');

U = Q * Uh(:,1:k);
S = S(1:k,1:k);
V = V(:,1:k);
end
