function [sol, sol1,SXX,SYY] = scca_init(X, Y, r, pen, tol, maxiter)

% Inputs:
% =======
% X, Y:    data sets, with rows correspond to samples and columns to variables
% r:       nuclear norm constraint
% pen:     penalty parameter on the l_1 norm of the solution, scaled by
%          sqrt(log(max(p1,p2))/n)
% tol:     tolerance level for convergence in ADMM
% maxiter: maximum number of iterations in ADMM
% 
% Outputs:
% ========
% sol:     optimum of the convex program
% sol1:    resacled optimum of the convex program, premultiplied by
%          SXX^(0.5), post multiplied by SYY^(0.5)

[n, ~] = size(X);

X = X - ones(n,1) * mean(X);
Y = Y - ones(n,1) * mean(Y);

SXX = X' * X / (n-1);
SXY = X' * Y / (n-1);
SYY = Y' * Y / (n-1);

[VX, DX] = eig(SXX);
DX = diag(DX);
idx = (abs(DX) > max(abs(DX)) * 1e-6);
DX = sqrt(DX(idx));
SXroot = VX(:,idx) * diag(DX) * (VX(:,idx)');
SXrootInv = VX(:,idx) * diag(1./DX) * (VX(:,idx)');

[VY, DY] = eig(SYY);
DY = diag(DY);
idy = (abs(DY) > max(abs(DY)) * 1e-6);
DY = sqrt(DY(idy));
SYroot = VY(:,idy) * diag(DY) * (VY(:,idy)');
SYrootInv = VY(:,idy) * diag(1./DY) * (VY(:,idy)');

% parameter in the augmented lagrangian
rho = 2;

[p1, p2] = size(SXY);

p = max(p1,p2);

B = SXrootInv * SXY * SYrootInv;


A = @(varargin)linop_crosprod( p1, p2, SXroot, SYroot, varargin{:} );

% x_cur = SXY ./ max(svd(SXY));
% x_cur = x_cur ./ max(sum(svd(x_cur)),r) .* r;

x_cur = 0;
y_cur = zeros(p1,p2);

[U D V] = svd(SXY, 'econ');
d = diag(D);
t = cap_soft_th(d, r, tol);
z_cur = U * diag(t) * V';

niter = 0;
initer = 0;

opts = [];
opts.printEvery = Inf;
opts.maxIts = 25;

while 1
    niter = niter + 1;
    initer = initer + 1;
    
    z_old = z_cur;
    Temp = x_cur - y_cur ./ rho +  B ./ rho;
    z_cur = tfocs(smooth_quad, {A, -Temp}, prox_l1(2 * pen * sqrt(log(p)/n) / rho), z_old, opts);
    
    x_old = x_cur;
    Temp = y_cur ./ rho + A(z_cur, 1);
    [U D V] = svd(Temp, 'econ');
    d = diag(D);
    t = cap_soft_th(d, r, tol);
    x_cur = U * diag(t) * V';
    
    y_old = y_cur;
    y_cur = y_old + rho .* (A(z_cur,1) - x_cur);
    
%     disp(niter);
    
    if max(rho*norm(x_cur - x_old, 'fro'), norm(z_cur - z_old, 'fro')) < tol
        break
    end
%     if initer == 100;
%         disp(max(norm(x_cur - A(z_cur,1)), rho.*norm(z_cur - z_old)));
%         initer = 0;
%     end
    if niter == maxiter
%         disp(max(norm(x_cur - A(z_cur,1)), rho.*norm(z_cur - z_old)));
        disp('Maximum number of iterations reached.');
        break
    end
end

sol = z_cur;
sol1 = x_cur;

% SCCALab v0.1 by Zongming Ma.
% All rights reserved.








