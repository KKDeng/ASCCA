function demo_random_cca(test_type)

if nargin<1; test_type = 1; end

isLasso = 0;


addpath ../SCCALab
addpath ../TFOCS-master
addpath ../
addpath ../util
addpath ../scca_penalty
addpath (genpath('../manopt'));

save_root = strcat('../results/random/');
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


N = 300; % sample size
p_set = [200,300]; % dim of X
q_set = [200,300]; % dim of Y
r_set = [1,2,3,4];   % column number
lambda_set = [0.1:0.1:0.8];
s  = [1:2:20];
ss = [1:4:20];
rng(1000); table_str = '';


for id_p =  1:length(p_set)        % n  dimension
    for id_q = 1 :length(q_set) % r  number of column
        for id_r = 1:length(r_set)  %mu  sparsity parameter
            bestTrace = 10; bestLasso = 10; bestInit = 10; bestPena = 10;
            for id_lambda = 1:length(lambda_set)
                lambda1 = lambda_set(id_lambda);
                lambda2 = lambda_set(id_lambda);
                r = r_set(id_r);
                p = p_set(id_p);
                q = q_set(id_q);
               
                
                
                basename = ['rand_',num2str(test_type),'_',num2str(p),'_',num2str(q),'_',num2str(r),'_',num2str(lambda1)];
                
                mu_X = zeros(p,1);   mu_Y = zeros(q,1);
                
                %% generate data matrix
                % Identity
                if test_type == 1
                    sigma_X = eye(p);   sigma_Y = eye(q);
                else
                    if test_type == 2
                        % toeplitz matrix
                        a = 0.3;
                        c1 = a.^((1:p)-1); c2 = a.^((1:q)-1);
                        sigma_X = toeplitz(c1);   sigma_Y = toeplitz(c2);
                    else
                        % sparse inverse
                        sigma_X = zeros(p); sigma_Y = zeros(q);
                        c1 =zeros(p,1); c1(1:3,:) = [1; 0.5; 0.4]; omegaX = toeplitz(c1);    sigma_X0 = inv(omegaX);
                        c2 =zeros(q,1); c2(1:3,:) = [1; 0.5; 0.4]; omegaY = toeplitz(c2);    sigma_Y0 = inv(omegaY);
                        for i=1:p
                            for j=1:p; sigma_X(i,j)= sigma_X0(i,j)/(sqrt(sigma_X0(i,i)*sigma_X0(j,j))); end
                        end
                        for i=1:q
                            for j=1:q; sigma_Y(i,j)= sigma_Y0(i,j)/(sqrt(sigma_Y0(i,i)*sigma_Y0(j,j))); end
                        end
                        sigma_X(abs(sigma_X)<1e-3)=0;  sigma_Y(abs(sigma_Y)<1e-3)=0;
                    end
                end
                
                if test_type == 4
                    sigma_X = eye(p);   sigma_Y = eye(q);
                    sigma_X(ss,ss) = 0.8; sigma_Y(ss,ss) = 0.8;
                    sigma_X = sigma_X - diag(diag(sigma_X)) + eye(p);
                    sigma_Y = sigma_Y - diag(diag(sigma_Y)) + eye(q);
                end
                
                if test_type == 5
                    sigma_X = eye(p);   sigma_Y = eye(q);
                    sigma_X(ss,ss) = 0.5; sigma_Y(ss,ss) = 0.5;
                    sigma_X = sigma_X - diag(diag(sigma_X)) + eye(p);
                    sigma_Y = sigma_Y - diag(diag(sigma_Y)) + eye(q);
                end
                
                if test_type == 6
                    sigma_X = eye(p);   sigma_Y = eye(q);
                    sigma_X(ss,ss) = 0.3; sigma_Y(ss,ss) = 0.3;
                    sigma_X = sigma_X - diag(diag(sigma_X)) + eye(p);
                    sigma_Y = sigma_Y - diag(diag(sigma_Y)) + eye(q);
                end
                % vector and rho
                u1= zeros(p,r);  v1 = zeros(q,r);
                
                
                lambda = diag( [0.9;0.8*ones(r-1,1);] );
                
                u1(s,(1:r)) = randi( [-2, 2], size(u1(s,(1:r))));
                Tss = sigma_X(s, s);
                u1 = u1 / sqrtm(u1(s,1:r)' * Tss * u1(s,1:r)+ 1e-13*eye(r));
                
                v1(s,(1:r)) = randi( [-2, 2], size(v1(s,(1:r))));
                Tss = sigma_Y(s, s);
                v1 = v1 / sqrtm(v1(s,1:r)' * Tss * v1(s,1:r)+ 1e-13*eye(r));
                
                
                [u_real, ~, ~] = svd(u1, 0);
                [v_real, ~, ~] = svd(v1, 0);
                %generate covariance matrix
                sigma_XY = sigma_X*u1*lambda*v1'*sigma_Y;
                
                
                
                
                Data =  mvnrnd([mu_X; mu_Y] ,[sigma_X,sigma_XY; sigma_XY', sigma_Y],10*N); % data matrix
                
                
                
                Xtrain = Data(1:N,1:p);     Xtest = Data(N+1:10*N,1:p);
                Ytrain = Data(1:N,p+1:p+q); Ytest = Data(N+1:10*N,p+1:p+q);
                
                
                
                
                [xhat,~,SXX,SYY] = scca_init(Xtrain, Ytrain, r,0.55, 1e-4, 10);
                [U0, ~, V0] = svds(xhat, r);
                Uinit_1 = U0(:,1:r);
                Vinit_1 = V0(:,1:r);
                [u0hat, ~,~] = svd(Uinit_1,'econ');  [v0hat,~,~] = svd(Vinit_1,'econ');
                [~,~,rho_init]  = canoncorr(Xtest * u0hat, Ytest * v0hat);
                lossu_init = norm(u0hat * u0hat'  - u_real * u_real', 'fro')^2;
                lossv_init = norm(v0hat * v0hat'  - v_real * v_real', 'fro')^2;
                
                
                %  CoLAR
                [Uhat, Vhat] = scca_refine(Xtrain,Ytrain, Uinit_1,Vinit_1, r, lambda1, lambda2);
                
                
                [X,Y,XtY,M1,M2] = normalize(Xtrain,Ytrain);
                
                [Up,Vp] = scca_penalty(M1,M2,XtY,X,Y,Uinit_1,Vinit_1,0.1*lambda1,0.1*lambda2);

                
                
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
                    h.cost = @(U,V,lambda) lambda*(sum(sum(abs(U))) + sum(sum(abs(V))));
                    h.prox = @proxLasso;
                    h.data = {0.1*lambda1};
                    
                end
                
                manifold = productmanifold(struct('U', stiefelgeneralizedfactory(p,r,M1),......
                    'V', stiefelgeneralizedfactory(q,r,M2)));
                
                % manifold = productmanifold(struct('U', stiefelfactory(p,r),......
                %     'V', stiefelfactory(q,r)));
                
                UV.U = Uinit_1/sqrtm(Uinit_1'*M1*Uinit_1);
                UV.V = Vinit_1/sqrtm(Vinit_1'*M2*Vinit_1);
                
                options_mialm.stepsize = 1/(2*abs(svds(full(XtY),1)));
                options_mialm.max_iter = 100;     options_mialm.maxitersub = 100;
                options_mialm.tau = 0.8;          options_mialm.rho = 1.05;
                options_mialm.nu0 = svds(X,1)^1*1 ; options_mialm.tol = 5e-4;
                options_mialm.gtol0 = 1;          options_mialm.gtol_decrease = 0.8;
                options_mialm.X0 = UV;      options_mialm.verbosity = 0;
                options_mialm.verbosity = 1;
                
                
                
                
                
                [X_mialm,Z1_mialm,Z2_mialm,out_mialm] = mialm(A,B, manifold, f, h, options_mialm);
                
                
                
                
                
                
                
                
                
                [u, ~, ~] = svd(X_mialm.U, 0);
                [v, ~, ~] = svd(X_mialm.V, 0);
                [Uhat, ~, ~] = svd(Uhat, 0);
                [Vhat, ~, ~] = svd(Vhat, 0);
                [Up, ~, ~] = svd(Up, 0);
                [Vp, ~, ~] = svd(Vp, 0);
                
                
                lossu_trace = norm(u * u'  - u_real * u_real', 'fro')^2;
                lossv_trace = norm(v * v'  - v_real * v_real', 'fro')^2;
                lossu_lasso = norm(Uhat * Uhat'  - u_real * u_real', 'fro')^2;
                lossv_lasso = norm(Vhat * Vhat'  - v_real * v_real', 'fro')^2;
                lossu_pena = norm(Up * Up'  - u_real * u_real', 'fro')^2;
                lossv_pena = norm(Vp * Vp'  - v_real * v_real', 'fro')^2;
                
                
                
                
                [~,~,rho_trace]  = canoncorr(Xtest * u, Ytest * v);
                [~,~,rho_lasso]  = canoncorr(Xtest * Uhat, Ytest * Vhat);
                [~,~,rho_pena]  = canoncorr(Xtest * Up, Ytest * Vp);
                
                
                save_path = strcat(save_root_res,basename,'.mat');
                save(save_path, 'lossu_trace', 'lossv_trace', 'lossu_pena', 'lossv_pena','lossu_lasso', 'lossv_lasso', 'lossu_init', 'lossv_init', ...
                    'rho_trace', 'rho_lasso', 'rho_init', 'rho_pena', 'Up','Vp', 'u', 'v', 'Uhat', 'Vhat', 'u0hat', 'v0hat','u_real','v_real', 'p','q','r', 'lambda1' );
                
                table_str = [table_str basename];
                table_str = [table_str sprintf('& %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f', ...
                    lossu_trace, lossv_trace,lossu_lasso,lossv_lasso, lossu_init, lossv_init,lossu_pena, lossv_pena)];
                table_str = [table_str '\\ \hline' newline];
                
                
                if bestTrace > mean([lossu_trace,lossv_trace])
                    bestTrace = mean([lossu_trace,lossv_trace]);
                    Trace = [lossu_trace,lossv_trace,rho_trace];
                    u_trace = u; v_trace = v;
                end
                if bestLasso > mean([lossu_lasso,lossv_lasso])
                    bestLasso = mean([lossu_lasso,lossv_lasso]);
                    Lasso = [lossu_lasso,lossv_lasso,rho_lasso];
                    u_lasso = Uhat; v_lasso = Vhat;
                end
                if bestInit > mean([lossu_init,lossv_init])
                    bestInit = mean([lossu_init,lossv_init]);
                    Init = [lossu_init,lossv_init,rho_init];
                end
                if bestPena > mean([lossu_pena,lossv_pena])
                    bestPena = mean([lossu_pena,lossv_pena]);
                    Pena = [lossu_pena,lossv_pena,rho_pena];
                    u_pena = Up;  v_pena = Vp;
                end
            end
            
            save_path = strcat(save_root_res1,basename,'.mat');
            save(save_path, 'Trace', 'Lasso', 'Init', 'Pena','p','q','r','u_trace','v_trace','u_lasso','v_lasso','u_pena','v_pena','u_real','v_real' );
        end
    end
end


disp(newline);
disp(table_str);
save_path = strcat(save_root,mat2str(test_type),'random.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',table_str);

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


