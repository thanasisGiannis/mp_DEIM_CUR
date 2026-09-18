
function [A] = load_mat_ex_1(m,n)
    A1 = zeros(m,n);
    A2 = zeros(m,n);
    for j = 1:10
        A1 = A1 + (2/j)*sprand(m,1,0.025)*sprand(n,1,0.025)';
    end
    for j = 11:300
        A2 = A2 + (1/j)*sprand(m,1,0.025)*sprand(n,1,0.025)';
    end
    A = A1+A2;
end