ss = RandStream('mt19937ar', 'Seed', 7);

addpath(genpath(cd));

disp('  ');
disp('--------------------------------------');
disp('Example: Toeplitz covariance matrix'); 
disp('n = 750, p1 = p2 = 200, s_u = s_v = 5');
disp('lambda1 = 0.9, lambda2 = 0.8');

n  = 750;
p1 = 200;
p2 = 200;
s  = [1, 6, 11, 16, 21];
lambda = diag( [0.9, 0.8] );
r = 2;

rr = 4;

disp('--------------------------------------');
disp('Generating data ...');

a = 0.3;
Sigma = eye(p1+p2);

T1 = toeplitz(a.^(0:1:(p1-1)));
Sigma(1:p1, 1:p1) = T1;
Tss = T1(s, s);
u = zeros(p1,r);
u(s,(1:r)) = randi(ss, [-2, 2], size(u(s,(1:r))));
u = u / sqrtm(u(s,1:r)' * Tss * u(s,1:r));


T2 = toeplitz(a.^(0:1:(p2-1)));
Sigma((p1+1):(p1+p2),(p1+1):(p1+p2)) = T2;
Tss = T2(s, s);
v = zeros(p2,r);
v(s,(1:r)) = randi(ss, [-2, 2], size(v(s,(1:r))));
v = v / sqrtm(v(s,1:r)' * Tss * v(s,1:r));


Sigma((p1+1):(p1+p2), 1:p1) = T2 *  v * lambda * u' * T1;
Sigma(1:p1, (p1+1):(p1+p2)) = Sigma((p1+1):(p1+p2), 1:p1)';

[u_n, ~, ~] = svd(u, 'econ');
[v_n, ~, ~] = svd(v, 'econ');

randn('state', 7);

Data = mvnrnd(zeros(p1+p2,1), Sigma, n);

X = Data(:,1:p1);
Y = Data(:,(p1+1):(p1+p2));

ntrain = floor(n/2);
Xtrain = X(1:ntrain, :);
Xtest = X((ntrain+1):n, :);
Ytrain = Y(1:ntrain, :);
Ytest = Y((ntrain+1):n, :);


disp('Data generated.');
disp('--------------------------------------');


tic;
% Initialization stage

[xhat, xhat_n] = scca_init(Xtrain, Ytrain, r, 0.55, 1e-4, 2e2); 

[U0, S0, V0] = svd(xhat, 'econ');
Uinit = U0(:,1:r);
Vinit = V0(:,1:r);

% Refinement stage
[uhat, vhat] = scca_refine(Xtrain, Ytrain, Uinit, Vinit, r, 1, 1);

[A B rho] = canoncorr(Xtest * uhat, Ytest * vhat);

disp('Canonical correlations on test data:');
disp(rho);

[uhat_n, ~, ~] = svd(uhat, 'econ');
[vhat_n, ~, ~] = svd(vhat, 'econ');
[Uinit, ~, ~] = svd(Uinit, 'econ');
[Vinit, ~, ~] = svd(Vinit, 'econ');

uerr = norm(uhat_n * uhat_n'  - u_n * u_n', 'fro');
verr = norm(vhat_n * vhat_n'  - v_n * v_n', 'fro');

uerr1 = norm(Uinit * Uinit'  - u_n * u_n', 'fro');
verr1 = norm(Vinit * Vinit'  - v_n * v_n', 'fro');

disp('Projection U error: ');
disp(uerr);
disp('Projection V error: ');
disp(verr);
toc;