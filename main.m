%% Optimization on the symplectic Stiefel manifold using BB-line searching
%
% We make use of the manifold structure in symplecticStiefelfactory
function main(M,S)
s = RandStream('mt19937ar'); % Random stream for reproducability

n = 100;
p = 5;
% E = [eye(p) zeros(p);
%             zeros(n-p,2*p);
%             zeros(p) eye(p);
%             zeros(n-p,2*p)];
% 
M = symplecticStiefelfactory(n,p,1);
problem.M = M;
% A = randn(s,2*n,2*p);
% A = A / norm(A,'fro');
% 
% 
% H = M.randham(n,s);
% H = H / norm(H,'fro');
%x0 = M.cay(H/2)*E;
% x0 = M.rand();
% f = @(x) 0.5*norm(x - A,'fro')^2;
% gf = @(x) (x - A);
% hf = @(x,v) 2*v;
U = randn(s,n,n) + randn(s,n,n)*i;
[U,~] = qr(U);

K = [real(U) -imag(U); imag(U) real(U)];

%A = 0.5 * (A + A');
% A = A*A';
% A = A / norm(A,2);
% % Initialize functions


l = n/5;
c = 1.2;
d = -sqrt(n/5);
L = gausstrans(n,l,c,d);


Q = K * L;
m = 2;
D = diag([zeros(m,1)' linspace(m+1,n,(n-m))]);

A = M.J(n) * Q * diag([diag(D)' diag(D)']) * (M.J(n) * Q )';
f = @(x) trace(x'*A * x);
gf = @(x) 2*A*x;
hf = @(x,v) 2*A*v;
% 
% r = 40;
% MAN = symplecticStiefelfactory(n,r,1);
% x = MAN.rand();
% v = MAN.randvec(x);
% [U,Sig,V] = svd(M.barOmega(x,v));
% 
% B = MAN.rand();
% c = rand(2*r,100);
% c = c / norm(c,'fro');
% S = B*c;
% 
% f = @(U) norm(S-U*M.Plus(U)*S,'fro')^2; % Objective function
% gf = @egradient; % The Euclidean gradient
% hf = @ehessian; % The Euclidean gradient

x0 = M.rand();

problem.M = M;
problem.cost = f;
problem.egrad = gf;
problem.ehess = hf;
options.gamma0 = f(x0);
options.gammamin = 1e-15;
options.gammamax = 1e15;
options.beta = 1e-4;
options.delta = 1e-1;
options.alpha = 0.85;
options.maxitter = 10000;
options.gtol = 1e-6;
options.mu = 10;
options.minstep = 1e-9;
tic;
[x1,info] = optim_sd(problem,x0,options);
toc; 
% problem.M = symplecticStiefelfactory(n,p,1);
% tic;
% [x2,info] = optim_sd(problem,x0,options);
% toc;
% Deig = eigs(A,M.J(n),2*n);
% tic;
options2.beta_type = 'S-D';
options2.verbosity = 2;
x = conjugategradient(problem,x0,options2);

%x = conjugategradient(problem,x0,options2);
x = trustregions(problem,x0);



function G = gausstrans(n,l,c,d)

    G = diag(ones(2*n,1));
    G(l-1,l-1) = c;
    G(l,l) = c;
    G(n+l-1,n+l-1) = 1/c;
    G(n+l,n+l) = 1/c;
    G(l,n+l-1) = d;
    G(l-1,n+l) = d; 
end
function ret = egradient(U)
    [n,~] = size(U);
    n = n/2;
    I = eye(2*n);
    ret = -2*((I-U*M.Plus(U))*S*S'*M.J(n)'*U*M.J(p) - M.J(n)*S*S'*(I-U*M.Plus(U))'*U*M.J(p));
end
function ret = ehessian(U,v)
    [n,~] = size(U);
    n = n/2;
    I = eye(2*n);
    SST = S * S';
    temp = -v*M.Plus(U);
    M1 = temp + M.Plus(temp);
    IUUPSSTJ = (I - U * M.Plus(U)) * SST * M.J(n)';

    
    M1SSTJT = M1 * SST * M.J(n)';
    %ret = -2 * ((M1SSTJT + M1SSTJT')* U * M.J(p) - (I - U * M.Plus(U)) *SST * M.J(n)'*v * M.J(p) - ((I - U * M.Plus(U)) *SST * M.J(n)')'*v*M.J(p)');
    
    ret = -2 * ((M1SSTJT- M1SSTJT')*U*M.J(p) + (IUUPSSTJ - IUUPSSTJ')* v * M.J(p)  );

end
end