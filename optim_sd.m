%% Steepest descend and conjugate gradient method using Riemannian BB line searching
function [x,info] = optim_sd(problem,x0,options)
tic;
M = problem.M;
x = x0;
xold = x;
info = cell(1,5);
gamma = options.gamma0;
c = problem.cost(x);
q = 1;
Z = -M.egrad2rgrad(x, problem.egrad(x));

for i = 1:options.maxitter
    %tic;
    
    

    if i > 1
        w = x - xold;
        y = Z - Zold;

        if mod(i-1,2) == 0
            gamma = norm(w,'fro')^2 / abs(trace(w'*y));
        else
            gamma = abs(trace(w'*y)) / norm(y,'fro')^2;
        end
    end
    
    gamma = max(options.gammamin,min(gamma,options.gammamax));
    
    j = 0;
    tau = gamma;
    
    INNER =  M.inner(x,-Z,Z);
    %toc;
    
    %tic;
    while true
        %tic;
        xnew = M.retr(x,Z,tau);
        cost = problem.cost(xnew);
        %toc;
        if cost <= c + options.beta * tau *INNER
            xold = x;
            x = xnew;
            Zold = Z;
            q = options.alpha * q + 1;
            c = (q - 1) / q * c + 1 / q * cost;
            break;
        else
            j = j + 1;
            tau = gamma * options.delta^j;
        end
        
    end
    Z = -M.egrad2rgrad(x, problem.egrad(x));
    gradnorm = M.norm(x,Z);
    info{i}.cost = cost;
    info{i}.gradnorm = gradnorm;
    info{i}.tau = tau;
    info{i}.iter = i;
    info{i}.time = toc;
    %toc;
    %tic;
    % fprintf(1, '%4d  %5.2e  %4.6e %3.2e\n', ...
    %      i, tau, cost, gradnorm);
        if gradnorm < options.gtol 
            fprintf('Algorithn converged!\n');
            info{i}.x = x;
            
            break;
        end
        if tau < options.minstep
            fprintf('Algorithm terminated due to small stepsize\n')
            info{i}.x = x;
            break;
        end
    %toc;
    

end

x = xnew;
end