function [Us,Ss,Vs] = my_svds_chop(A,k,tol,copt)


A = chop(A,copt);
[m,n] = size(A);

S = 0;
u_prev = chop(zeros(m,1),copt);

beta   = 0;
betas  = [];
alphas = [];

V = [];
U = [];
%AV = [];
B  = [];

v = chop(randn(n,1),copt);
v = chop(v./chop(norm(v),copt),copt);
for jj=1:min(m,n)
    for j=max(size(V,2),1): max(3*k,15)
        u = chop(chop(A*v,copt)-chop(beta*u_prev,copt),copt);
        if j>1
            u = chop(u-chop(U*(chop(U'*u,copt)),copt),copt);
            u = chop(u-chop(U*(chop(U'*u,copt)),copt),copt);
        end
        alpha = chop(norm(u),copt);
        alphas = [alphas alpha];

        u = chop(u/alpha,copt);
        V = [V v];
        %AV = [AV A*v];
        v = chop(chop(A'*u,copt)-chop(alpha*v,copt),copt);

        if j>1
            v = chop(v-chop(V*(chop(V'*v,copt)),copt),copt);
            v = chop(v-chop(V*(chop(V'*v,copt)),copt),copt);
        end

        beta = chop(norm(v),copt);
        betas  = [betas beta];
        v = chop(v./beta,copt);

        U = [U u];
        u_prev = u;
        %B = [B U(:,1:end-1)'*AV(:,end); U(:,end)'*AV(:,1:end-1) U(:,end)'*AV(:,end)];
        B = diag(alphas) + diag(betas(1:end-1),1);
        if j>=k
            %[Us,Ss,Vs] = svd(U'*AV,'econ');
            [P,Ss,Q] = svd(B,'econ');

            P = chop(P(:,1:k),copt);
            Q = chop(Q(:,1:k),copt);
            Ss = chop(Ss(1:k,1:k),copt);

            Us = chop(U*P,copt);
            Vs = chop(V*Q,copt);

            nrmResEst = beta * max(abs(P(end,1:k)));
            if nrmResEst < tol
                Ss = diag(Ss);
                return;
            end
        end
    end
    
    V = Vs;
    U = Us;
    AV = chop(A*V,copt);
    B  = chop(U'*AV,copt);
    betas  = [0 diag(B,1)']; 
    alphas = diag(B)';
end
Ss = diag(Ss);
end
