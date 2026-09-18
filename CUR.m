
function [C,U,R,eta_p,eta_q] = CUR(A,V,W,options)
    [J] = DEIM(W, options); 
    C = A(:,J);
    
    [I] = DEIM(V,options);
    R = A(I,:);
    
    U = C\A/R;
    
    eta_q = norm(inv(W(J,:)));
    eta_p = norm(inv(V(I,:)));

end

function p = DEIM(V,options)

    chop([],options)
    
    V = chop(V);
    [m,k]=size(V);
    
    v = chop(V(:,1));
    [~,p1]=max(abs(v));
    p = p1;
    for j = 2:k
        v = chop(V(:, j));
        c = chop(chop(V(p, 1:j-1)) \ chop(v(p)));
        r = chop(chop(v) - chop(V(:,1:j-1))*chop(c));
        [~, pj] = max(abs(r));
        p = [p, pj];
    end
end