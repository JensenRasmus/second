%% Numerical examples 
function numerical_examples(M,S)
clear all;
s = RandStream('mt19937ar');
RandStream.setGlobalStream(s);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                          Numerical Experiments                          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Experiment 1: Nearest symplectic matrix problem
% Experiment 2: Symplectic eigenvalue problem (Sterming from Bin-Gao 2024)
% Experiment 3: Proper orthogonal decompoosition of low-rank matrix

%% Experiment 1: Nearest symplectic matrix problem

n = 1000;
ks = [10,50,100];

num_iter_r_sd = [];
runtime_r_sd = [];
grad_norm_r_sd = [];
feasibility_r_sd = [];
f_val_r_sd = [];

num_iter_r_cg = [];
runtime_r_cg = [];
grad_norm_r_cg = [];
feasibility_r_cg = [];
f_val_r_cg = []; 

num_iter_r_tr = [];
runtime_r_tr = [];
grad_norm_r_tr = [];
feasibility_r_tr = [];
f_val_r_tr = [];

for k = 1:3

    p = ks(k);

    M = symplecticStiefelfactory(n,p,1);
    A = rand(2*n,2*p);
    A = A / norm(A,'fro');
    f = @(x) 1/2 * norm(A - x,'fro')^2;
    gf = @(x) x-A;
    hf = @(x,v) v;
    problem.M = M;
    problem.cost = f;
    problem.egrad = gf;
    problem.ehess = hf;

    x0 = M.rand();

    options.gamma0 = f(x0);
    options.gammamin = 1e-15;
    options.gammamax = 1e15;
    options.beta = 1e-4;
    options.delta = 1e-1;
    options.alpha = 0.85;
    options.maxitter = 10000;
    options.gtol = 1e-6;
    options.mu = 5;
    options.minstep = 1e-11;


    [x_sd,info_sd] = optim_sd(problem,x0,options);
    [x_cg,info_cg] = optim_cg(problem,x0,options);
    [x_tr,cost,info_tr] = trustregions(problem,x0);

    runtime_r_sd(k) = info_sd{end}.time;
    num_iter_r_sd(k) = info_sd{end}.iter;
    grad_norm_r_sd(k) = info_sd{end}.gradnorm;
    feasibility_r_sd(k) = M.checkmanifold(info_sd{end}.x);
    f_val_r_sd(k) = f(x_sd);

    runtime_r_cg(k) = info_cg{end}.time;
    num_iter_r_cg(k) = info_cg{end}.iter;
    grad_norm_r_cg(k) = info_cg{end}.gradnorm;
    feasibility_r_cg(k) = M.checkmanifold(info_cg{end}.x);
    f_val_r_cg(k) = f(x_cg);

    runtime_tr = [info_tr.time];
    runtime_r_tr(k) = runtime_tr(end);

    num_iter_tr = [info_tr.iter];
    num_iter_r_tr(k) = num_iter_tr(end);

    grad_norm_tr = [info_tr.gradnorm];
    grad_norm_r_tr(k) = grad_norm_tr(end);

    feasibility_r_tr(k) = M.checkmanifold(x_tr);
    f_val_r_tr(k) = f(x_tr);
end

varNames = ["itterations", "runtime (s)", "Grad norm at end", "Feasibility", "f(x^*)"];
rowNames = ["10","50","100"];
T1 = table(num_iter_r_sd',runtime_r_sd',grad_norm_r_sd',feasibility_r_sd',f_val_r_sd', 'rownames',rowNames,'Variablenames',varNames);

T2 = table(num_iter_r_cg',runtime_r_cg',grad_norm_r_cg',feasibility_r_cg',f_val_r_cg', 'rownames',rowNames,'Variablenames',varNames);

T3 = table(num_iter_r_tr',runtime_r_tr',grad_norm_r_tr',feasibility_r_tr',f_val_r_tr', 'rownames',rowNames,'Variablenames',varNames);

disp("########")
disp("# R-SD #")
disp("########")
disp(T1);

disp("########")
disp("# R-CG #")
disp("########")
disp(T2);

disp("########")
disp("# R-TR #")
disp("########")
disp(T3);

%% Experiemnt 2: The symplectic eigenvalue problem
% n = 100;
% p = 5;
% M = symplecticStiefelfactory(n,p,1);
% 
% rng default
% U = randn(s,n,n) + randn(s,n,n)*i;
% %U = randn(n,n) + randn(n,n)*i;
% [U,~] = qr(U); % Unitary matrix U are obtained from QR
% 
% K = [real(U) -imag(U); imag(U) real(U)];
% 
% l = n/5;
% c = 1.2;
% d = -sqrt(n/5);
% L = gausstrans(n,l,c,d);
% 
% Q = K * L;
% m = 2;
% %D = diag([zeros(m,1)' linspace(m+1,n,(n-m))]);
% D = diag(linspace(1,n,n));
% A = M.J(n) * Q * diag([diag(D)' diag(D)']) * (M.J(n) * Q )';
% %A = Q * diag([diag(D)' diag(D)'] )*Q';
% f = @(x) trace(x'*A * x);
% gf = @(x) 2*A*x;
% hf = @(x,v) 2*A*v;
% 
% problem.M = M;
% problem.cost = f;
% problem.egrad = gf;
% problem.ehess = hf;
% 
% x0 = M.rand();
% 
% options.gamma0 = f(x0);
% options.gammamin = 1e-15;
% options.gammamax = 1e15;
% options.beta = 1e-4;
% options.delta = 1e-1;
% options.alpha = 0.85;
% options.maxitter = 10000;
% options.gtol = 1e-6;
% options.mu = 5;
% options.minstep = 1e-11;
% 
% [x_sd,info_sd] = optim_sd(problem,x0,options);
% [x_cg,info_cg] = optim_cg(problem,x0,options);
% [x_tr,cost,info_tr] = trustregions(problem,x0);
% 
% 
% runtime_r_sd = info_sd{end}.time;
% num_iter_r_sd = info_sd{end}.iter;
% grad_norm_r_sd = info_sd{end}.gradnorm;
% feasibility_r_sd = M.checkmanifold(info_sd{end}.x);
% error_sd = abs(f(x_sd) - 30);
% 
% runtime_r_cg = info_cg{end}.time;
% num_iter_r_cg = info_cg{end}.iter;
% grad_norm_r_cg = info_cg{end}.gradnorm;
% feasibility_r_cg = M.checkmanifold(info_cg{end}.x);
% error_cg = abs(f(x_cg) - 30);
% 
% runtime_tr = [info_tr.time];
% runtime_r_tr = runtime_tr(end);
% 
% num_iter_tr = [info_tr.iter];
% num_iter_r_tr = num_iter_tr(end);
% 
% grad_norm_tr = [info_tr.gradnorm];
% grad_norm_r_tr = grad_norm_tr(end);
% 
% feasibility_r_tr = M.checkmanifold(x_tr);
% 
% error_tr = abs(f(x_tr) - 30);
% 
% [S,D_sd,Q,G] = williasondiag(x_sd'*A*x_sd);
% [S,D_cg,Q,G] = williasondiag(x_cg'*A*x_cg);
% [S,D_tr,Q,G] = williasondiag(x_tr'*A*x_tr);
% 
% varNames = ["itterations", "runtime (s)", "Grad norm at end", "Feasibility", "||f(x^*)-30||"];
% T1 = table(num_iter_r_sd',runtime_r_sd',grad_norm_r_sd',feasibility_r_sd',error_sd', 'Variablenames',varNames);
% T2 = table(num_iter_r_cg',runtime_r_cg',grad_norm_r_cg',feasibility_r_cg',error_cg', 'Variablenames',varNames);
% T3 = table(num_iter_r_tr',runtime_r_tr',grad_norm_r_tr',feasibility_r_tr',error_tr', 'Variablenames',varNames);
% 
% disp("########")
% disp("# R-SD #")
% disp("########")
% %disp(T1);
% 
% disp("########")
% disp("# R-CG #")
% disp("########")
% disp(T2);
% 
% disp("########")
% disp("# R-TR #")
% disp("########")
% disp(T3);
% format long
% disp("Computed eigenvalues")
% disp("   R-SD                R-CG                R_TR")
% disp(horzcat(D_sd,D_cg,D_tr))


%% Experiment 3: The proper symplectic decomposition
% % 
% n = 1000;
% ks = [10,20,40];
% 
% num_iter_r_sd = [];
% runtime_r_sd = [];
% grad_norm_r_sd = [];
% feasibility_r_sd = [];
% f_val_r_sd = [];
% 
% num_iter_r_cg = [];
% runtime_r_cg = [];
% grad_norm_r_cg = [];
% feasibility_r_cg = [];
% f_val_r_cg = [];
% 
% num_iter_r_tr = [];
% runtime_r_tr = [];
% grad_norm_r_tr = [];
% feasibility_r_tr = [];
% f_val_r_tr = [];
% 
% for k = 1:3
% 
%     p = ks(k);
%     M = symplecticStiefelfactory(n,p,1);
% 
%     r = 40;
%     MAN = symplecticStiefelfactory(n,r,1);
% 
%     B = MAN.rand();
%     c = rand(2*r,100);
%     c = c / norm(c,'fro');
%     S = B*c;
%     f = @(x) norm(S - x*M.Plus(x)*S,'fro')^2;
%     gf = @egradient;
%     hf = @ehessian;
%     problem.M = M;
%     problem.cost = f;
%     problem.egrad = gf;
%     problem.ehess = hf;
% 
%     x0 = M.rand();
% 
%     options.gamma0 = f(x0);
%     options.gammamin = 1e-15;
%     options.gammamax = 1e15;
%     options.beta = 1e-4;
%     options.delta = 1e-1;
%     options.alpha = 0.85;
%     options.maxitter = 10000;
%     options.gtol = 1e-6;
%     options.mu = 5;
%     options.minstep = 1e-11;
% 
%     [x_sd,info_sd] = optim_sd(problem,x0,options);
%     [x_cg,info_cg] = optim_cg(problem,x0,options);
%     [x_tr,cost,info_tr] = trustregions(problem,x0);
% 
%     runtime_r_sd(k) = info_sd{end}.time;
%     num_iter_r_sd(k) = info_sd{end}.iter;
%     grad_norm_r_sd(k) = info_sd{end}.gradnorm;
%     feasibility_r_sd(k) = M.checkmanifold(info_sd{end}.x);
%     f_val_r_sd(k) = f(x_sd);
% 
%     runtime_r_cg(k) = info_cg{end}.time;
%     num_iter_r_cg(k) = info_cg{end}.iter;
%     grad_norm_r_cg(k) = info_cg{end}.gradnorm;
%     feasibility_r_cg(k) = M.checkmanifold(info_cg{end}.x);
%     f_val_r_cg(k) = f(x_cg);
% 
%     runtime_tr = [info_tr.time];
%     runtime_r_tr(k) = runtime_tr(end);
% 
%     num_iter_tr = [info_tr.iter];
%     num_iter_r_tr(k) = num_iter_tr(end);
% 
%     grad_norm_tr = [info_tr.gradnorm];
%     grad_norm_r_tr(k) = grad_norm_tr(end);
% 
%     feasibility_r_tr(k) = M.checkmanifold(x_tr);
%     f_val_r_tr(k) = f(x_tr);
% end
% 
% varNames = ["itterations", "runtime (s)", "Grad norm at end", "Feasibility", "||S - U^*(U^*)^+S||"];
% rowNames = ["10","20","40"];
% T1 = table(num_iter_r_sd',runtime_r_sd',grad_norm_r_sd',feasibility_r_sd',f_val_r_sd', 'rownames',rowNames,'Variablenames',varNames);
% T2 = table(num_iter_r_cg',runtime_r_cg',grad_norm_r_cg',feasibility_r_cg',f_val_r_cg', 'rownames',rowNames,'Variablenames',varNames);
% T3 = table(num_iter_r_tr',runtime_r_tr',grad_norm_r_tr',feasibility_r_tr',f_val_r_tr', 'rownames',rowNames,'Variablenames',varNames);
% 
% disp("########")
% disp("# R-SD #")
% disp("########")
% disp(T1);
% 
% disp("########")
% disp("# R-CG #")
% disp("########")
% disp(T2);
% 
% disp("########")
% disp("# R-TR #")
% disp("########")
% disp(T3);


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
    
    ret = -2 * ((M1SSTJT- M1SSTJT')*U*M.J(p) + (IUUPSSTJ - IUUPSSTJ')* v * M.J(p)  );

end
end