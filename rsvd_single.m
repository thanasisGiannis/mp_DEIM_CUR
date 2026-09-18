function [U,S,V] = rsvd_single(A,k,p,q)
%RSVD_SINGLE  Rank-k randomized SVD in single precision.
%
% All explicitly created arrays and matrix products use single precision.
% The precision actually used by QR/SVD is MATLAB-release and platform
% dependent; check with class(U), class(S), and class(V).

if nargin < 3 || isempty(p), p = 10; end
if nargin < 4 || isempty(q), q = 0;  end

A = single(A);
[m,n] = size(A);
ell = min(k + p, min(m,n));

% Ensure that the random test matrix is single.
Omega = randn(n, ell, 'single');
Y = A * Omega;

% Optional subspace iteration.
for j = 1:q
    [Q,~] = qr(Y, 0);
    Y = A' * Q;
    [Q,~] = qr(Y, 0);
    Y = A * Q;
end

[Q,~] = qr(Y, 0);

B = Q' * A;
[Uh,S,V] = svd(B, 'econ');

U = Q * Uh(:,1:k);
S = S(1:k,1:k);
V = V(:,1:k);
end
