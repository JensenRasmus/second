function Om = barOmega(U,Delta)
    [n,~] = size(U);
    J = @(n) [zeros(n,n), eye(n); -eye(n), zeros(n,n)];
    
    M = U/(U'*U);
    X = J(n/2)'*M*Delta'*J(n/2)';
    Om = Delta*M' + X - X*M*U';
end