function M = symplecticStiefelfactory(n,p,type)
% Returns a manifold struct to optimize over the symplectic Stiefel
% manifold SpSt(2n,2p). Default: gpuflag = 'false', type = 1
% Call as
% 
% M = symplecticStiefelfactory(n,p)
% M = symplecticStiefelfactory(n,p,type)
%
% to obtain M in SpSt(2n,2p).
% 
% The metric is the right-invariant metric obtained from 
% Bendokat & Zimmermann 2021.
% 
% Points are represented by matrices U of size 2n x 2p 
% 
% The points have to fulfill the conditon that Plus(U)*U = eye(2p), where
% Plus(U) = J(2p)' * U' * J(2n), where 
%
% J(2n) = |  0    eye(n)|
%         |-eye(n   0   |.
%
% Tangent vectors are reperesented as matrices of size 2n x 2p.
%
% Type: 
% type == 0: The retraction is taken to approximate the geodesics wrt. the 
%           pseudo-Riemannian metric in Bendokat & Zimmermann 2021. 
%           The vector transport is taken as the differentiated retraction. 
% 
% type == 1: The retraction is taken to approximate the geodesics wrt. the 
%           right-invariant metric in Bendokat & Zimmermann.
%           The vector transport is taken as the differentiated retraction. 
%
% type == 2: The retraction is taken to approximate the geodesics wrt. the 
%           right-invariant metric in Bendokat & Zimmermann. The vector
%           transport is taken as the orthogonal projection.
%
% Author: Rasmus Jensen


    % Checks and initialization of standard config
    assert(n >= p, 'The dimension n must be larger than the dimension p.');
    
    if ~exist('type','var') || isempty(type)
        type = 1;
    end


    % Not used due to numerical instability. 
    % Enable to use, for example, 'checkhessian'
    M.exp = @geod;
    function geo = geod(U,eta,t)
        %eta = eta / norm(eta,'fro');
        omega = M.barOmega(U,eta);
        geo = expm(t*(omega-omega'))*expm(t*omega')*U;
    end

    M.name = @() sprintf('symplectic Stiefel manifold SpSt(%d,%d)',2*n,2*p);
    
    M.dim = @() (4*n-2*p+1) * p;
    
    M.J = @(d) [zeros(d) eye(d);
              -eye(d) zeros(d)];

    M.Plus = @Plus;
    function PlusA = Plus(A)
        % Calculate the symplectic inverse A^+ = -JA'J of A
        % where A = [a1, a2; a3, a4]
        %     J = @(n) [zeros(n,n), eye(n); -eye(n), zeros(n,n)];
        [s,t] = size(A);
        if mod(s,2) == 0 && mod(t,2) == 0
            a1 = A(1:s/2,1:t/2);
            a2 = A(1:s/2,t/2+1:t);
            a3 = A(s/2+1:s,1:t/2);
            a4 = A(s/2+1:s,t/2+1:t);
            PlusA = [a4', -a2'; -a3', a1'];
    %        AP = -J(k/2)*A'*J(n/2);
        else
            error("Plus requires a 2n-by-2k matrix as input.")
        end
    end

    M.checktangent = @(U,eta) norm(M.Plus(eta)*U +M.Plus(U)*eta,'fro');
    
    M.inner = @inp;
    
    function res = inp(U,eta1,eta2)
        UTU = (U'*U);
        UTUinv = inv(U'*U);
        res = trace( (eta1'*( eta2 - 0.5*M.Plus(U')* (M.Plus(UTUinv)*(M.Plus(U)*eta2)) ))*UTUinv );
    end
    %M.norm = @(U,eta) norm(eta,'fro'); % Using the Frobenious norm
    M.norm = @(U,eta) sqrt(M.inner(U,eta,eta));

    M.dist = @(x, y) error('stiefel.dist not implemented yet.');

    M.typicaldist = @() sqrt((4*n-2*p+1) * p);

    M.tangent2ambient_is_identity = true;
    M.tangent2ambient = @(X, U) U;
    
    % Convert Euclidean gradients to Riemannian gradients
    
    M.egrad2rgrad = @egrad2rgrad;
    function rgrad = egrad2rgrad(U,egrad)
        %rgrad = egradUT * U + M.J(n) * egradUT' * M.J(n) * U;
        rgrad = egrad*(U'*U) - M.Plus(U') * (M.Plus(egrad) * U);
    end
    M.barOmega = @barOmega;
    function Om = barOmega(U,Delta)
        Q = U/(U'*U);
        
        X = -M.Plus(Delta*Q');
        Om = Delta*Q' + X - (X*Q)*U';

    end
    %M.ehess2rhess = @ehesstorhess;
    M.ehess2rhess = @ehesstorhess2;
    function g = Gamma(U,eta)
    % Compute the Christoffel function when the twp inputs are the same
        % tic;
        omega = M.barOmega(U,eta);
        % toc;
        %g = -(omega-omega') * (eta + omega'*U)-(omega')^2*U;
        OTU = omega' * U;
        
        g = -(omega-omega') * (eta + OTU)-omega'*OTU;
        
        % tic;
        % g2 = (-omega + omega' ) * eta - omega * (OTU);
        % toc
        % norm(g - g2)
    end
    function hess = ehesstorhess(U,eg,eh,v)
        % For manopt routines this should be enabled. They are essentially
        % just identity maps when eg and eh are inputted as matrices.
        eg = @(U) eg;
        eh = @(U,v) eh;

        egrad = eg(U);
        ehess = eh(U,v);
        
        grad = M.egrad2rgrad(U,egrad);
        if norm(grad - v) < eps
           norm_gradmv = 1; 
        else
            norm_gradmv = norm(grad - v,'fro');
        end

        N1 = norm((grad +v),'fro');
        
        % FYY = ehess * (U' * U) + egrad*v' * U + (egrad)* U' * v + ...
        %       M.J(n) * v * (egrad)' * M.J(n) * U + ...
        %       M.J(n) * U *(ehess)' * M.J(n) * U + ...
        %       M.J(n) * U * (egrad)' * M.J(n) * v;
        
        VTU = v' * U; % Update May 24'
        
        FYY = ehess * (U'*U) + egrad*VTU + (egrad)* VTU' + ...
            -(M.Plus(egrad*v')  + M.Plus(ehess*U'))*U - M.Plus(egrad * U') * v;
        
        hess = FYY + 1/4 * (N1^2*Gamma(U,(grad +v)/N1)-norm_gradmv^2*Gamma(U,(grad-v)/norm_gradmv));
    end
    
    function approxHess = ehesstorhess2(U,eg,eh,v)
        eg = @(U) eg;
        eh = @(U,v) eh;
        
        egrad = eg(U);
        ehess = eh(U,v);
 
        VTU = v' * U; 
        
        FYY = ehess * (U'*U) + egrad*VTU + (egrad)* VTU' + ...
            -(M.Plus(egrad*v')  + M.Plus(ehess*U'))*U - M.Plus(egrad * U') * v;
    
        approxHess = M.tangent(U,FYY);
    end
    M.retr_right_inv = @retraction_right_inv;
    function [ret,storage] = retraction_right_inv(U,eta,t)
        if ~exist('t', 'var')|| isempty(t)
            t = 1;
        end
        
        UTUinv = eye(2*p) / (U' *U);
        UPlus = M.Plus(U);
        UUPlus =  U * UPlus;

        A = (-M.Plus(U) * M.Plus(eta')) * M.Plus(UTUinv)+UTUinv * (eta' * U) - UTUinv * (eta' *M.Plus(U') * M.Plus(UTUinv))  ;
        H = -(eye(2*n) - UUPlus)*(M.Plus(eta')*M.Plus(UTUinv));

        barDelta = U * A + H; % This is the horizontal lift of eta
        PlusbarDelta = M.Plus(barDelta);
        X = [((eye(2 * n) - 0.5 * UUPlus) * barDelta) -U];
        Y = [(UPlus') (PlusbarDelta * (eye(2*n) - 0.5 * UUPlus))'];
        
        hatX = [Y -X];
        hatY = [X Y];
        hatYTX = hatY' * hatX; 
        w_to_save = t* hatX / (eye(8*p) - t/2 * hatYTX) * hatY';
        w = eye(2*n) + w_to_save;
        theta = (eye(2 * p) - t/2 * A + t^2/4 * M.Plus(H) * H);
        cay_right = (-U + (t*H + 2*U) / theta);
        ret = w * cay_right;
        if nargout == 2
            % Store for computing a vector transport
            storage.o1hat = w_to_save;
            storage.t = t;
            storage.X = X;
            storage.Y = Y;
            storage.hatX = hatX;
            storage.hatY = hatY;
            storage.hatYTX = hatYTX;
            storage.cay_right = cay_right;
        end
    end
    
    M.retr_pseudo_rie = @retraction_pseudo_rie;
    function ret = retraction_pseudo_rie(U,eta,t)
        if ~exist('t','var') || isempty(t)
            t = 1;
        end
        A = M.Plus(U)*eta;
        H = eta - U*A;
        I = t^2/4*M.Plus(H)*H - t/2*A + eye(2*p);
        ret = -U + (t*H+2*U) * inv(I);
        
    end
    
    % Vector transport associated with retr_right_inv
    
    M.transp_right_inv = @trans_right_inv;
    function res = trans_right_inv(U,xi,teta,storage)
        if ~exist('storage','var') || isempty(storage)
            t = 1;
            % No previous information can be aquired
            %storage = 0;
            % Compute without cashing
            UTUinv = eye(2*p) / (U' *U);
            UPlus = M.Plus(U);
            UUPlus =  U * UPlus;
            teta = teta / t;

            A = (-M.Plus(U) * M.Plus(teta')) * M.Plus(UTUinv)+UTUinv * (teta' * U) - UTUinv * (teta' *M.Plus(U') * M.Plus(UTUinv))  ;
            H = -(eye(2*n) - UUPlus)*(M.Plus(teta')*M.Plus(UTUinv));
            
            barDelta = U * A + H; % This is the "pseudo-Riemannaian tangent representation" of teta
            PlusbarDelta = M.Plus(barDelta);
            X_teta = [((eye(2 * n) - 0.5 * UUPlus) * barDelta) -U];
            Y_teta = [(UPlus') (PlusbarDelta * (eye(2*n) - 0.5 * UUPlus))'];
            
            hatX_teta = [Y_teta -X_teta];
            hatY_teta = [X_teta Y_teta];
            hatYTX_teta = hatY_teta' * hatX_teta; 

            %omega_teta = Y_teta*X_teta';

            % A2 = M.J(p) * U' * xi * UTUinv * M.J(p) + UTUinv * (xi' * U) - UTUinv * xi' * (M.J(n)' * U * UTUinv * M.J(p));
            % H2 = (eye(2*n) - UUPlus) * M.J(n) * xi * UTUinv * M.J(p);
    
            A2 = (-M.Plus(U) * M.Plus(xi')) * M.Plus(UTUinv)+UTUinv * (xi' * U) - UTUinv * (xi' *M.Plus(U') * M.Plus(UTUinv))  ;
            H2 = -(eye(2*n) - UUPlus)*(M.Plus(xi')*M.Plus(UTUinv));
            

            barDelta = U * A2 + H2; % This is the "pseudo-Riemannaian tangent representation" of xi
            PlusbarDelta = M.Plus(barDelta);
            X_xi = [((eye(2 * n) - 0.5 * UUPlus) * barDelta) -U];
            Y_xi = [(UPlus') (PlusbarDelta * (eye(2*n) - 0.5 * UUPlus))'];
            hatX_xi = [Y_xi -X_xi];
            hatY_xi = [X_xi Y_xi];
            
            O1_xi = hatX_xi * hatY_xi' ;

            omega_xi = Y_xi*X_xi';
            cay_right_part = -U + (t*H + 2*U) / (eye(2*p) - 0.5 * t * A + t^2 * 1/4 * M.Plus(H)*H);

            M1_temp = 0.5 * t * hatX_teta / (eye(8*p) - 0.5 * t * hatYTX_teta) * hatY_teta';
            M1 = eye(2*n) + M1_temp; 
            M2 = cay_right_part;
            M3 = 2*M1_temp + eye(2*n); 
            M4 = eye(2*n) + 0.5 *t* X_teta / (eye(4*p) - 0.5 *t* Y_teta' * X_teta) * Y_teta';

            res_temp = (M1 * 0.5 * (O1_xi * (M1 * M2))) + ...
                       (M3 * (M4 * 0.5*(omega_xi' * (M4 * U))));
            res = norm(xi,'fro')/norm(res_temp,'fro') * res_temp;

        else 
            assert(class(storage) == "struct", "The last input has to either be empty or a storage structure")

            w_to_save = storage.o1hat;
            t = storage.t;
            X = storage.X;
            Y = storage.Y;
            hatX = storage.hatX;
            hatY = storage.hatY;
            hatYTX = storage.hatYTX;
            cay_right = storage.cay_right;

            M12 = eye(2*n) + 0.5 * w_to_save;
            O1_xi2 = hatX * hatY';
            
            P1 = (M12 * 0.5 * (O1_xi2 * (M12 * cay_right)));
            M22 = eye(2*n) + w_to_save;
            O2_xi2 = Y*X';
            M32 = eye(2*n) + 0.5 * t * X / (eye(4*p) - t/2*Y'*X)*Y';

            P2 = (M22 * (M32 * 0.5 * (O2_xi2' * (M32 * U))));

            res_temp = P1 + P2;
            res = norm(xi,'fro') / norm(res_temp,'fro') * res_temp;
            % norm(M1 - M12)
            % norm(O1_xi - O1_xi2)
            % norm(P1 - M1 * 0.5 * O1_xi * M1 * M2)
            % norm(M3 - M22)
            % norm(omega_xi - O2_xi2)
            % norm(M32 - M4)
            % norm(P2 - M3 * M4 * 0.5*omega_xi' * M4 * U)
            % 
            % norm((P1 + P2) - res_temp)

            s = 0;
            % t = 1;
            % % No previous information can be aquired
            % %storage = 0;
            % % Compute without cashing
            % UTUinv = eye(2*p) / (U' *U);
            % UPlus = M.Plus(U);
            % UUPlus =  U * UPlus;
            % teta = teta / t;
            % 
            % A = (-M.Plus(U) * M.Plus(teta')) * M.Plus(UTUinv)+UTUinv * (teta' * U) - UTUinv * (teta' *M.Plus(U') * M.Plus(UTUinv))  ;
            % H = -(eye(2*n) - UUPlus)*(M.Plus(teta')*M.Plus(UTUinv));
            % 
            % barDelta = U * A + H; % This is the "pseudo-Riemannaian tangent representation" of teta
            % PlusbarDelta = M.Plus(barDelta);
            % X_teta = [((eye(2 * n) - 0.5 * UUPlus) * barDelta) -U];
            % Y_teta = [(UPlus') (PlusbarDelta * (eye(2*n) - 0.5 * UUPlus))'];
            % 
            % hatX_teta = [Y_teta -X_teta];
            % hatY_teta = [X_teta Y_teta];
            % hatYTX_teta = hatY_teta' * hatX_teta; 
            % 
            % %omega_teta = Y_teta*X_teta';
            % 
            % % A2 = M.J(p) * U' * xi * UTUinv * M.J(p) + UTUinv * (xi' * U) - UTUinv * xi' * (M.J(n)' * U * UTUinv * M.J(p));
            % % H2 = (eye(2*n) - UUPlus) * M.J(n) * xi * UTUinv * M.J(p);
            % 
            % A2 = (-M.Plus(U) * M.Plus(xi')) * M.Plus(UTUinv)+UTUinv * (xi' * U) - UTUinv * (xi' *M.Plus(U') * M.Plus(UTUinv))  ;
            % H2 = -(eye(2*n) - UUPlus)*(M.Plus(xi')*M.Plus(UTUinv));
            % 
            % 
            % barDelta = U * A2 + H2; % This is the "pseudo-Riemannaian tangent representation" of xi
            % PlusbarDelta = M.Plus(barDelta);
            % X_xi = [((eye(2 * n) - 0.5 * UUPlus) * barDelta) -U];
            % Y_xi = [(UPlus') (PlusbarDelta * (eye(2*n) - 0.5 * UUPlus))'];
            % hatX_xi = [Y_xi -X_xi];
            % hatY_xi = [X_xi Y_xi];
            % 
            % O1_xi = hatX_xi * hatY_xi' ;
            % 
            % omega_xi = Y_xi*X_xi';
            % cay_right_part = -U + (t*H + 2*U) / (eye(2*p) - 0.5 * t * A + t^2 * 1/4 * M.Plus(H)*H);


            % In the storage structure
            % X_xi = storage.X
            % Y_xi = storage.Y
            % hatX_xi = hatX_teta (when teta = teta / t) above 
            % hatY_xi = hatY_teta (when teta = teta / t) above
            % t = storage.t
            % teta = t*xi
            % M2 = storage.cay_right

            % Forming the matrices from the storage structure:
            % hatX_teta = [storage.Y -storage.X];
            % hatY_teta = [storage.X storage.Y];
            % hatYTX_teta = hatY_teta' * hatX_teta; 
            % t = storage.t;
            % X_teta = storage.X;
            % Y_teta = storage.Y;
            % cay_right_part = storage.cay_right;
            % 
            % O1_xi = hatX_teta * hatY_teta';
            % omega_xi = Y_teta * X_teta';
            
            % hatX_teta2 = [storage.Y -storage.X];
            % hatY_teta2 = [storage.X storage.Y];
            % hatYTX_teta2 = hatY_teta2' * hatX_teta2; 
            % t2 = storage.t;
            % X_teta2 = storage.X;
            % Y_teta2 = storage.Y;
            % cay_right_part2 = storage.cay_right;
            % 
            % O1_xi2 = hatX_teta2 * hatY_teta2';
            % omega_xi2 = Y_teta2 * X_teta2';
            % 
        end
            
        

        % M1_temp2 = 0.5 * t2 * hatX_teta2 / (eye(8*p) - 0.5 * t2 * hatYTX_teta2) * hatY_teta2';
        % M12 = eye(2*n) + M1_temp2; 
        % M22 = cay_right_part2;
        % M32 = 2*M1_temp2 + eye(2*n); 
        % M42 = eye(2*n) + 0.5 *t2* X_teta2 / (eye(4*p) - 0.5 *t2* Y_teta2' * X_teta2) * Y_teta2';

        
        % res_temp2 = M12 * 0.5 * O1_xi2 * M12 * M22 + ...
        %                M32 * M42 * 0.5*omega_xi2' * M42 * U;
        % M1 * 0.5 * O1_xi * M1 * M2  and M12 * 0.5 * O1_xi2 * M12 * M22
        % differ in a strange way. 
        % res2 = norm(xi,'fro')/norm(res_temp2,'fro') * res_temp2;
        % norm(res - res2,'fro')

    end
    
    M.transp_pseudo_rie = @trans_pseudo_rie;
    function res = trans_pseudo_rie(U,xi,teta,storage)
        if ~exist('storage','var') || isempty(storage)
            storage = 0;
        end
        OmegatEta = M.barOmega(U,teta);
        OmegaXi = M.barOmega(U,xi);
        I = eye(2*n);

        res = inv(I - 0.5*OmegatEta) * 0.5 * OmegaXi * inv(I - 0.5 * OmegatEta) * U;
    end
    switch type
        case 0
            % Retraction wrt. the pseudo-Riemannian metric
            M.retr = M.retr_pseudo_rie;
            M.transp = M.transp_pseudo_rie;
        case 1
            % Retraction wrt. the right-invariant metric
            % (Default choice)
            M.retr = M.retr_right_inv;
            %M.retr = M.proj;
            M.transp = M.transp_right_inv;
        case 2
            M.retr = M.retr_right_inv;
            M.transp = @proj;
    end
    

    % This is taken form stiefelfactory.m, and the author have
    % not considered its relevance
    M.hash = @(X) ['z' hashmd5(X(:))];
    

    % Generating random matrices is done via the Calyey transform
    % Cay:hamil(2n) -> Sp(2n), 
    % where hamil(2n) is the space of Hamiltonian matrices, and Sp(2n) is
    % the Symplectic group. Mapping from Sp(2n) to the symplectic Stiefel
    % manifold SpSt(2n,2k) is done via the quotiemt map pi(M)=ME, where
    % E = | eye(k)  0    |
    %     |   0     0    |
    %     |   0   eye(k) |
    %     |   0     0    |  
    %
    % It can be useful to generate random Hamiltonian matrices, and therefore
    % the the function M.randhamil(n) is supplied.
    %
    % We also supply functions to compute the Cayley retraction and 
    % horizontal lift of a tangent.  

    M.randham = @randomhamil;
    function Omega = randomhamil(n,s)
        if nargin < 2 || isempty(s)
            A = randn(n);
            B = randn(n);
            C = randn(n);
        else
            A = randn(s,n);
            B = randn(s,n);
            C = randn(s,n);
        end
        
        B = (B + B')/2;
        C = (C + C')/2;
        
        Omega = [A, B; C, -A'];
    end
    
    M.cay = @cayley;
    function cay = cayley(X)
        cay = (eye(2*n) + X) / (eye(2*n) - X);
    end


    M.rand = @randompoint;
    function point = randompoint()
        Omega = M.randham(n);
        Omega = Omega / norm(Omega,'fro');
        E = [eye(p) zeros(p);
            zeros(n-p,2*p);
            zeros(p) eye(p);
            zeros(n-p,2*p)];
        point = M.cay(Omega/2)*E;
    end
    
    % Generate random tangent vector
    % M.randvec = M.egrad2rgrad;
    M.randvec = @generatetangent;
    function tan = generatetangent(U)
        Omega = M.randham(n);
        tan = Omega * U;
        %tan / norm(tan,'fro');
    end
    
    % M.tangent = @project;
    % function ret = project(x,u)
    %     % Make sure that u is a tangnet vector by projecting 
    %     % it onto the tangent space.
    %     % Taken from the source code of the Julia implementation of manopt
    %     % @online{2106.08777,
    %     %     Author = {Seth D. Axen and Mateusz Baran and Ronny Bergmann and Krzysztof Rzecki},
    %     %     Title = {Manifolds.jl: An Extensible Julia Framework for Data Analysis on Manifolds},
    %     %     Year = {2021},
    %     %     Eprint = {2106.08777},
    %     %     Eprinttype = {arXiv},
    %     % }
    %     Q = M.J(n) * x;
    % 
    %     uTQ = u' * Q;
    %     W = uTQ - uTQ';
    %     xTx = x' * x;
    %     L = sylvester(xTx,xTx,-0.5 * W);
    %     ret = u - Q * (L - L');
    % end

    M.tangent = @proj;
    function ret = proj(U,v,a1,a2)
        % The fields a1 and a2 are just placeholders, ensuring a flexible
        % implementation.
        
        if ~isfield('a1','var') || isempty(a1)
            a1 = 0;
        end
        if ~isfield('a2','var') || isempty(a2)
            a2 = 0;
        end
        % U is the point on the manifold, v is the vector which is
        % projected.
        if n > p
            UTP = M.Plus(U');
            UTUPinv = M.Plus(U'*U);
            vPlusU = M.Plus(v)*U;

            ret = v - 0.5 * UTP * (UTUPinv \ (vPlusU + M.Plus(vPlusU )) );
        else
            vPlusU = M.Plus(v)*U;

            ret = v - 0.5 * U * (vPlusU + M.Plus(vPlusU));
        end
    end
    
    M.checkmanifold = @checkmanifold;
    function res = checkmanifold(U)
        res = norm(M.Plus(U) * U - eye(2*p));
    end

    M.lincomb = @matrixlincomb;

    M.zerovec = @(x) zeros(2*n, 2*p);
end