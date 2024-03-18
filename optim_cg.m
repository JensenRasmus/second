%% Steepest descend and conjugate gradient method using Riemannian BB line searching
function [x,info] = optim_cg(problem,x0,options)
tic;
M = problem.M;
x = x0;
xold = x;
info = cell(1,5);

gamma = options.gamma0;
c = problem.cost(x);
q = 1;
pk = -M.egrad2rgrad(x, problem.egrad(x));
gradk = pk;
gradkp1 = -gradk;

for i = 1:options.maxitter
    %tic;
    
    if i > 1
        w = x - xold;
        y = gradkp1 - gradk;

        if mod(i-1,2) == 0
            gamma = norm(w,'fro')^2 / abs(trace(w'*y));
        else
            gamma = abs(trace(w'*y)) / norm(y,'fro')^2;
        end
    end
    
    gamma = max(options.gammamin,min(gamma,options.gammamax));
    
    j = 0;
    tau = gamma;
    
    INNER =  M.inner(x,-gradk,pk);
    %toc;
    
    %tic;
    while true
        %tic;
        [xnew,storage] = M.retr(x,pk,tau);
        % xnew = M.retr(x,pk,tau);
        
        cost = problem.cost(xnew);
        %toc;
        if cost <= c + options.beta * tau *INNER
            xold = x;
            x = xnew;
            %Zold = Z;
            q = options.alpha * q + 1;
            c = (q - 1) / q * c + 1 / q * cost;
            break;
        else
            j = j + 1;
            tau = gamma * options.delta^j;
        end
        
    end
    gradk = gradkp1;
    gradkp1 = M.egrad2rgrad(x,problem.egrad(x));
    

    gradnorm = M.norm(x,gradkp1);
    if mod(i,options.mu) == 0
        pk = -gradkp1;
    else
        num = gradnorm^2;
        %l = 1;
        %num = gradnorm^2 - M.inner(x,gradkp1,l * M.transp(xold,gradk,tau*pk));
        denom = M.inner(xold,gradk,gradk);
        betak = num/denom;

        tpk = M.transp(xold,pk,tau*pk,storage);

        % tpk = M.transp(xold,pk,tau*pk);
        % norm(tpk - tpk2)
        pk = -gradkp1 + betak * tpk;

    end

    
    info{i}.cost = cost;
    info{i}.gradnorm = gradnorm;
    info{i}.tau = tau;
    info{i}.iter = i;
    info{i}.time = toc;
    %toc;
    %tic;
    % fprintf(1, '%4d  %5.2e  %4.6e %3.2e\n', ...
    %     i, tau, cost, gradnorm);
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


end