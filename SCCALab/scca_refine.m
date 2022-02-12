function [Uhat, Vhat] = scca_refine(X, Y, U0, V0, r, penU, penV)

% Inputs:
% =======
% X, Y:      data sets, with rows correspond to samples and columns to variables
% U0, V0:    initial estimators for U and V
% r:         desired rank
% penU,penV: penalty parameter on the l_2-l_1 norm of the estimators, scaled by
%            sqrt( (r + log(max(p1,p2))) / n )
% 
% Outputs:
% ========
% Uhat,Vhat: estimators for U and V
% 

opts = [];
opts.printEvery = Inf;


[n, p1] = size(X);
[~, p2] = size(Y);

X = X - ones(n,1) * mean(X);
Y = Y - ones(n,1) * mean(Y);

p = max(p1,p2);

X = X./sqrt(n);
Y = Y./sqrt(n);

Uhat = tfocs(smooth_quad, {linop_matrix(X, 'R2R', r), - Y*V0}, prox_l1l2(penU*sqrt((r+log(p))/n)), U0, opts);
Vhat = tfocs(smooth_quad, {linop_matrix(Y, 'R2R', r), - X*U0}, prox_l1l2(penV*sqrt((r+log(p))/n)), V0, opts);


% a = 1;


% SCCALab v0.1 by Zongming Ma.
% All rights reserved.
