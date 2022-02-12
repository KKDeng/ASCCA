function demo_realdata(test_type)

if nargin<1; test_type = 2; end

isLasso = 0;


addpath ../SCCALab
addpath ../TFOCS-master
addpath ../
addpath ../util
addpath ../scca_penalty
addpath ../data
addpath (genpath('../manopt'));

save_root = strcat('../results/real/');
if ~exist(save_root,'dir')
    mkdir(save_root);
end


save_root_res = strcat(save_root,'res/');
if ~exist(save_root_res,'dir')
    mkdir(save_root_res);
end

save_root_res1 = strcat(save_root,'res_s/');
if ~exist(save_root_res1,'dir')
    mkdir(save_root_res1);
end



r = test_type;
lambda_set = [0.1:0.1:0.8];

rng(1000); table_str = '';


data = load('geneData.mat');
for id_num = 21:23  %mu  sparsity parameter
    bestTrace = 0; bestLasso = 0; bestInit = 0; bestPena = 0;
    DataX = data.A{id_num,1}'; DataY = data.B{id_num,1}';
    [N,p] = size(DataX); q = size(DataY,2); n = ceil(N/4)*3;
    Xtrain = DataX(1:n,:);     Xtest = DataX(n+1:N,:);
    Ytrain = DataY(1:n,:);     Ytest = DataY(n+1:N,:);
    N = n;
    [X,Y,XtY,M1,M2] = normalize(Xtrain,Ytrain);
    for id_lambda = 1:length(lambda_set)
        lambda1 = lambda_set(id_lambda);
        lambda2 = lambda1;
         basename = ['rand_',num2str(test_type),'_',num2str(p),'_',num2str(q),'_',num2str(r),'_',num2str(lambda1)];
    
        [xhat,~,SXX,SYY] = scca_init(Xtrain, Ytrain, r,0.55, 1e-4, 10);
        [U0, ~, V0] = svds(xhat, r);
        Uinit_1 = U0(:,1:r);
        Vinit_1 = V0(:,1:r);
        [u0hat, ~,~] = svd(Uinit_1,'econ');  [v0hat,~,~] = svd(Vinit_1,'econ');
        [~,~,rho_init]  = canoncorr(Xtest * u0hat, Ytest * v0hat);
       
        
        
        %  CoLAR
        [Uhat, Vhat] = scca_refine(Xtrain,Ytrain, Uinit_1,Vinit_1, r, lambda1, lambda2);
        
        
        
        A = struct();
        A.applyA = @AXu;
        A.applyAT = @AtXu;
        
        B = struct();
        B.applyB = @BYv;
        B.applyBT = @BtYu;
        
        f = struct();
        f.cost_grad = @cca_cost_grad;
        f.data = {XtY};
        
        h = struct();
        h.cost = @(U,V,lambda1,lambda2) lambda1*norm(svd(U),1) + lambda2*norm(svd(V),1);
        h.prox = @proxNuclear;
        h.data = {0.15*lambda1,0.15*lambda2};
        
        
        if isLasso==1
            A = struct();
            A.applyA = @(X) X;
            A.applyAT = @(X) X;
            
            B = struct();
            B.applyB = @(X) X;
            B.applyBT = @(X) X;
            h = struct();
            h.cost = @(U,V,lambda1,lambda2) lambda1*sum(sum(abs(U))) + lambda2*sum(sum(abs(V)));
            h.prox = @proxLasso;
            h.data = {0.15*lambda1,0.15*lambda2};
            
        end
        
        manifold = productmanifold(struct('U', stiefelgeneralizedfactory(p,r,M1),......
            'V', stiefelgeneralizedfactory(q,r,M2)));
        
        % manifold = productmanifold(struct('U', stiefelfactory(p,r),......
        %     'V', stiefelfactory(q,r)));
        
        UV.U = Uinit_1/sqrtm(Uinit_1'*M1*Uinit_1);
        UV.V = Vinit_1/sqrtm(Vinit_1'*M2*Vinit_1);
        
        options_mialm.stepsize = 1/(2*abs(svds(full(XtY),1)));
        options_mialm.max_iter = 30;     options_mialm.maxitersub = 100;
        options_mialm.tau = 0.8;          options_mialm.rho = 1.05;
        options_mialm.nu0 = svds(X,1)^1*1 ; options_mialm.tol = 1e-3;
        options_mialm.gtol0 = 1;          options_mialm.gtol_decrease = 0.8;
        options_mialm.X0 = UV;      options_mialm.verbosity = 0;
        options_mialm.verbosity = 1;
        
        
        
        
        
        [X_mialm,Z1_mialm,Z2_mialm,out_mialm] = mialm(A,B, manifold, f, h, options_mialm);
        
        
        [Up,Vp] = scca_penalty(M1,M2,XtY,X,Y,Uinit_1,Vinit_1,0.1*lambda1,0.1*lambda2);
        
        
        
        
        
        
        [u, ~, ~] = svd(X_mialm.U, 0);
        [v, ~, ~] = svd(X_mialm.V, 0);
        [Uhat, ~, ~] = svd(Uhat, 0);
        [Vhat, ~, ~] = svd(Vhat, 0);
        [Up, ~, ~] = svd(Up, 0);
        [Vp, ~, ~] = svd(Vp, 0);
        
        
       
        
        
        
       
        
        
        [~,~,rho_trace]  = canoncorr(Xtest * u, Ytest * v);
        [~,~,rho_lasso]  = canoncorr(Xtest * Uhat, Ytest * Vhat);
        [~,~,rho_pena]  = canoncorr(Xtest * Up, Ytest * Vp);
        
        
        save_path = strcat(save_root_res,basename,'.mat');
        save(save_path, 'rho_trace', 'rho_lasso', 'rho_pena', 'rho_init', 'Up', 'Vp', 'u', 'v', 'Uhat', 'Vhat', 'u0hat', 'v0hat', 'p','q','r', 'lambda1' );
        
        
        
        
        if bestTrace < mean(rho_trace)
            bestTrace = mean(rho_trace);
            Trace = rho_trace;
            lambda_trace = lambda1;
            u_trace = u; v_trace = v;
            [~,~,rho_trace_t]  = canoncorr(Xtrain * u, Ytrain * v);
        end
        if bestLasso < mean(rho_lasso)
            bestLasso = mean(rho_lasso);
            Lasso = rho_lasso;
            lambda_lasso = lambda1;
            u_lasso = Uhat; v_lasso = Vhat;
            [~,~,rho_lasso_t]  = canoncorr(Xtrain * Uhat, Ytrain * Vhat);
        end
        if bestInit < mean(rho_init)
            bestInit = mean(rho_init);
            Init = rho_init;
            lambda_init = lambda1;
            u_init = u0hat; v_init = v0hat;
             [~,~,rho_init_t]  = canoncorr(Xtrain * u0hat, Ytrain * v0hat);
        end
        if bestPena < mean(rho_pena)
            bestPena = mean(rho_pena);
            Pena = rho_pena;
            lambda_pena = lambda1;
            u_pena = Up; v_pena = Vp;
             [~,~,rho_pena_t]  = canoncorr(Xtrain * Up, Ytrain * Vp);
        end
    end
    
    save_path = strcat(save_root_res1,basename,'.mat');
    save(save_path, 'Trace', 'Lasso', 'Init','Pena', 'p','q','r',...
        'lambda_trace','u_trace','v_trace','rho_trace_t','lambda_lasso','u_lasso','v_lasso','rho_lasso_t',...
        'lambda_pena','u_pena','v_pena','rho_pena_t','lambda_init','u_init','v_init','rho_init_t');
end



% disp(newline);
% disp(table_str);
% save_path = strcat(save_root,'random.txt');
% fid = fopen(save_path,'w+');
% fprintf(fid,'%s',table_str);

    function AX = AXu(U)
        AX = zeros(N,p*r);
        for k=1:r
            AX(:,(k-1)*p+1:k*p) = X.*repmat(U(:,k)',N,1);
        end
    end


    function BY = BYv(V)
        BY = zeros(N,q*r);
        for k=1:r
            BY(:,(k-1)*q+1:k*q) = Y.*repmat(V(:,k)',N,1);
        end
    end

    function AtX = AtXu(W)
        AtX = zeros(p,r);
        for k=1:r
            AtX(:,k) = sum(X.*W(:,(k-1)*p+1:k*p))';
        end
    end


    function BtY = BtYu(W)
        BtY = zeros(q,r);
        for k=1:r
            BtY(:,k) = sum(Y.*W(:,(k-1)*q+1:k*q))';
        end
    end



    function [f,g] = cca_cost_grad(UV,XtY)
        g.U = -XtY*UV.V;  g.V = -XtY'*UV.U;
        f = sum(sum((g.V).*(UV.V)));
    end



    function g = proxNuclear(U,V,mu,lambda1,lambda2)
        [g.U,~] = prox_nuclear(U,lambda1*mu);
        [g.V,~] = prox_nuclear(V,lambda2*mu);
    end

    function g = proxLasso(U,V,mu,lambda1,lambda2)
        g.U = max(abs(U) - mu*lambda1,0).* sign(U);
        
        g.V = max(abs(V) - mu*lambda2,0).* sign(V);
    end



end


